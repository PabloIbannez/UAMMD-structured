#pragma once

#include "misc/IBM.cuh"
#include "misc/IBM_kernels.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct ExternalTabulated_{

        using Kernel = uammd::IBM_kernels::Peskin::threePoint;

        //Storage data
        struct StorageData{

                thrust::device_vector<real3> forceData;
                thrust::device_vector<real> energyData;

                std::shared_ptr<thrust::device_vector<real>> energy_i;
                std::shared_ptr<thrust::device_vector<real3>> force_i;

                std::shared_ptr<uammd::IBM<Kernel>> interpolator;
        };

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   const int& nx,const int& ny,const int& nz,
                                                   const real& scale,
                                                   const std::vector<int>& id_x,
                                                   const std::vector<int>& id_y,
                                                   const std::vector<int>& id_z,
                                                   const std::vector<real>& energy,
                                                   const std::vector<real3>& force) {

                Box box = gd->getEnsemble()->getBox();

                System::log<System::MESSAGE>("[External tabulated] Number of cells: (%d,%d,%d).", nx, ny, nz);

                System::log<System::MESSAGE>("[External tabulated] Reading external tabulated potential.");
                if (nx <= 0 || ny <= 0 || nz <= 0) {
                    System::log<System::CRITICAL>("[External tabulated] Cell dimensions must be positive.");
                }

                Grid grid(box, make_int3(nx, ny, nz));

                System::log<System::MESSAGE>("[External tabulated] Grid cell size: (%f,%f,%f).",
                                              grid.cellSize.x, grid.cellSize.y, grid.cellSize.z);

                real h = grid.cellSize.x;
                if (std::abs(grid.cellSize.x-grid.cellSize.y) > 1e-6 ||
                    std::abs(grid.cellSize.x-grid.cellSize.z) > 1e-6) {
                    System::log<System::CRITICAL>("[External tabulated] Grid cell size must be isotropic.");
                } else {
                    System::log<System::MESSAGE>("[External tabulated] Grid cell size is isotropic. Using h = %f.", h);
                }

                auto kernel = std::make_shared<Kernel>(h);


                if (id_x.size() != nx*ny*nz) {
                    System::log<System::CRITICAL>("[External tabulated] The number of given data values does not match the number of cells.");
                }


                // Rearrange energy vector, this is done to match the ibm layout
                thrust::host_vector<real>  rearranged_energy(nx*ny*nz, real(0.0f)); // Initialize with zeros
                thrust::host_vector<real3> rearranged_force( nx*ny*nz, make_real3(0.0f, 0.0f, 0.0f)); // Initialize with zeros
                for (size_t i = 0; i < id_x.size(); ++i) {
                    UAMMD_SET_3D_ROW_MAJOR(rearranged_energy, nz, ny, nx, id_z[i], id_y[i], id_x[i], scale*energy[i]);
                    UAMMD_SET_3D_ROW_MAJOR(rearranged_force,  nz, ny, nx, id_z[i], id_y[i], id_x[i], scale*force[i]);
                }

                //Set up storage data

                StorageData storage;

                storage.energyData = rearranged_energy;
                storage.forceData  = rearranged_force;

                storage.energy_i = std::make_shared<thrust::device_vector<real>>(pg->getNumberParticles());
                thrust::fill(storage.energy_i->begin(),storage.energy_i->end(),0.0);

                storage.force_i = std::make_shared<thrust::device_vector<real3>>(pg->getNumberParticles());
                thrust::fill(storage.force_i->begin(),storage.force_i->end(),make_real3(0.0));

                storage.interpolator = std::make_shared<uammd::IBM<Kernel>>(kernel,grid);

                return storage;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data) {


                int nx = data.getParameter<int>("nx");
                int ny = data.getParameter<int>("ny");
                int nz = data.getParameter<int>("nz");

                real scale = 1.0;
                if (data.isParameterAdded("scale")) {
                    scale = data.getParameter<real>("scale");
                    System::log<System::MESSAGE>("[External tabulated] Scaling factor (energy and force): %f.", scale);
                } else {
                    System::log<System::MESSAGE>("[External tabulated] No scaling factor provided."
                                                 " Assuming energy and force are in the correct units.");
                }

                std::vector<int> id_x  = data.getData<int>("i");
                std::vector<int> id_y  = data.getData<int>("j");
                std::vector<int> id_z  = data.getData<int>("k");

                std::vector<real> energy = data.getData<real>("energy");
                std::vector<real3> force = data.getData<real3>("force");

                return getStorageData(gd,pg,nx,ny,nz,scale,id_x,id_y,id_z,energy,force);
        };

        //Computational data
        struct ComputationalData{
            Box box;

            real4* pos;

            //Result vectors
            real* energy_i;
            real3* force_i;
        };


        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData& storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st) {

            Box box = gd->getEnsemble()->getBox();

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            ComputationalData computational;

            computational.box      = box;
            computational.pos      = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.energy_i = thrust::raw_pointer_cast(storage.energy_i->data());
            computational.force_i  = thrust::raw_pointer_cast(storage.force_i->data());

            if (comp.energy == true or comp.lambdaDerivative == true) {
                storage.interpolator->gather(computational.pos,
                                             computational.energy_i,
                                             thrust::raw_pointer_cast(storage.energyData.data()),
                                             pd->getNumParticles(),st);
                //cudaStreamSynchronize(st);
                cudaDeviceSynchronize();
            }

            if (comp.force == true) {
                storage.interpolator->gather(computational.pos,
                                             computational.force_i,
                                             thrust::raw_pointer_cast(storage.forceData.data()),
                                             pd->getNumParticles(),st);
                //cudaStreamSynchronize(st);
                cudaDeviceSynchronize();
            }

            return computational;
        }


        static inline __device__ real energy(int index_i,const ComputationalData& computational){
                real e = computational.energy_i[index_i];
                computational.energy_i[index_i] = real(0.0);
                return e;
        }

        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
                real3 f = computational.force_i[index_i];
                computational.force_i[index_i] = make_real3(0.0);
                return f;
        }

    };

    using ExternalTabulated   = External_<ExternalTabulated_>;

}}}}

