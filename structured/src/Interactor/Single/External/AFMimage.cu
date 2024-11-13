#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/External/External.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Definitions/Matrices.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct AFMimage_{

        //Storage data
        struct StorageData{
            thrust::device_vector<real> Hexp;
        };


        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data) {

            real dx_img = data.getParameter<real>("dx"); //Pixel size in x
            real dy_img = data.getParameter<real>("dy"); //Pixel size in y

            std::vector<int> px  = data.getData<int>("x"); // Pixel x coordinate
            std::vector<int> py  = data.getData<int>("y"); // Pixel y coordinate

            size_t Nx = px.size(); // Number of pixels in x direction
            size_t Ny = py.size(); // Number of pixels in y direction

            std::vector<real> height = data.getData<real>("height"); // Height of the pixel

            /////////////////////////////
            // INPUT DATA VERIFICATION //
            /////////////////////////////

            // Get box and compare the pixel size with the given values of dx and dy
            Box box = gd->getEnsemble()->getBox();

            real dx = box.boxSize.x / Nx;
            real dy = box.boxSize.y / Ny;

            if(real(abs(dx - dx_img)) > real(1e-6) || real(abs(dy - dy_img)) > real(1e-6)){
                System::log<System::CRITICAL>("[AFMimage] The pixel size (%f,%f) does not match the given values of dx and dy (%f,%f)",
                                              dx,dy,dx_img,dy_img);
            }

            // Check if the pixel coordinates are within the box
            for(uint i = 0; i < Nx; i++){
                if(px[i] < 0 || px[i] >= Nx){
                    System::log<System::CRITICAL>("[AFMimage] Detected pixel x coordinate outside the box: %d (should be in [0,%d])",
                                                  px[i],Nx);
                }
            }

            for(uint i = 0; i < Ny; i++){
                if(py[i] < 0 || py[i] >= Ny){
                    System::log<System::CRITICAL>("[AFMimage] Detected pixel y coordinate outside the box: %d (should be in [0,%d])",
                                                  py[i],Ny);
                }
            }

            // Check all pixels are given. For 0,...,Nx-1 and 0,...,Ny-1 all combinations should be present
            std::vector<bool> pixels(Nx*Ny,false);
            for(uint i = 0; i < height.size(); i++){
                int px_i = px[i];
                int py_i = py[i];

                UAMMD_SET_2D_ROW_MAJOR(pixels,Ny,Nx,py_i,px_i,true);
            }

            for(uint i = 0; i < Nx; i++){
                for(uint j = 0; j < Ny; j++){
                    bool present = UAMMD_GET_2D_ROW_MAJOR(pixels,Ny,Nx,j,i);
                    if(!present){
                        System::log<System::CRITICAL>("[AFMimage] Pixel (%d,%d) is missing",i,j);
                    }
                }
            }

            // Check if the height is positive
            for(uint i = 0; i < height.size(); i++){
                if(height[i] < 0){
                    System::log<System::CRITICAL>("[AFMimage] Detected negative height value: %f (should be positive)",
                                                  height[i]);
                }
            }

            //////////////////
            // STORAGE DATA //
            //////////////////

            StorageData storage;

            thrust::host_vector<real> Hexp_h(Nx*Ny,0.0);

            for(uint i = 0; i < height.size(); i++){
                int px_i = px[i];
                int py_i = py[i];

                // #define UAMMD_SET_2D_ROW_MAJOR(matrix, rows, cols, row, col, value) (matrix)[(row) * (cols) + (col)] = (value)
                // The number of rows is Ny and the number of columns is Nx
                UAMMD_SET_2D_ROW_MAJOR(Hexp_h.data(),Ny,Nx,py_i,px_i,height[i]);
            }

            storage.Hexp = Hexp_h; // Data is copied to the GPU

        };

        //Computational data
        struct ComputationalData{
            Box box;

            real4* pos;
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

            return computational;
        }


        static inline __device__ real energy(int index_i,const ComputationalData& computational){
            real e = real(0.0);
            return e;
        }

        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
            real3 f = make_real3(0.0);
            return f;
        }

    };

    using AFMimage   = External_<AFMimage_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    External,AFMimage,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::AFMimage>
)

