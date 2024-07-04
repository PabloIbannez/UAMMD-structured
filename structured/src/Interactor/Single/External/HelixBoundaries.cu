#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/External/External.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <cub/cub.cuh>

#include "ExternalTabulated.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    namespace HelixBoundaries_ns{

        inline constexpr uint threadsPerBlock = 256;

        inline __device__ real3 helixEquation(real t, real pitch, real radius, real eps,Box box){
            real x = radius*cos(t);
            real y = eps*radius*sin(t); // if eps = 1, the helix is right-handed, if eps = -1, the helix is left-handed
            real z = pitch*t/real(2*M_PI);
            z = z - box.boxSize.z/real(2.0);
            return make_real3(x,y,z);
        }

        inline __device__ real constraintPotential(real r, real R, real K){
                // r is the distance from a point to the helix axis
                // R is the helix inner radius
                if(r < R){
                   return real(0.0);
                }
                return real(0.5)*K*(r - R)*(r - R);
        }

        template<int threadsPerBlock>
        __global__ void fillGridDistanceEnergy2Helix(real helixPitch,
                                                     real helixRadius,
                                                     real eps,
                                                     real helixInnerRadius,
                                                     real K,
                                                     real tMax,
                                                     Box box,
                                                     real h, // cell size
                                                     int blockIterationsPerPoint,
                                                     int nPointsHelix,
                                                     real* __restrict__ energy) {
            // Calculate current cell index
            const uint3 cellIndex = blockIdx;
            real3 cellPos = make_real3(cellIndex.x * h - box.boxSize.x / real(2.0) + h / real(2.0),
                                       cellIndex.y * h - box.boxSize.y / real(2.0) + h / real(2.0),
                                       cellIndex.z * h - box.boxSize.z / real(2.0) + h / real(2.0)); // Center of the cell

            typedef cub::BlockReduce<real, threadsPerBlock> BlockReduce;
            __shared__ typename BlockReduce::TempStorage temp_storage;

            real threadData;

            // Initialize minDistance with a large value
            real minDistance = box.boxSize.x;

            for (int i = 0; i < blockIterationsPerPoint; i++) {
                int pointIndex = i * blockDim.x + threadIdx.x;
                real t = real(pointIndex) / real(nPointsHelix); // t in [0, 1]
                t = t * tMax; // t in [0, tMax]

                if (pointIndex < nPointsHelix) {
                    real3 dr = box.apply_pbc(helixEquation(t, helixPitch, helixRadius, eps, box) - cellPos);
                    threadData = sqrt(dot(dr, dr));
                } else {
                    // We want to discard these values
                    threadData = box.boxSize.x; // Since we are looking for the minimum.
                }

                // Perform block-wide reduction to find the minimum distance
                real blockMin = BlockReduce(temp_storage).Reduce(threadData, cub::Min());

                if (threadIdx.x == 0) {
                    minDistance = min(minDistance, blockMin);
                }
                __syncthreads();
            }

            // Now thread 0 has the minimum distance for this block
            if (threadIdx.x == 0) {
                real pot = constraintPotential(minDistance, helixInnerRadius, K);
                UAMMD_SET_3D_ROW_MAJOR(energy, gridDim.z, gridDim.y, gridDim.x, cellIndex.z, cellIndex.y, cellIndex.x, pot);
            }
        }

        __global__ void computeForce(Box    box,
                                     Grid   grid,
                                     real*  energy,
                                     real3* gradient) {
            // One thread per cell
            int idx = blockIdx.x*blockDim.x + threadIdx.x;
            int nCells = grid.getNumberCells();

            if(idx < nCells) {

                real h = grid.cellSize.x; // Assuming cubic cells

                int nX = grid.cellDim.x;
                int nY = grid.cellDim.y;

                int i = idx % nX;
                int j = (idx / nX) % nY;
                int k = idx / (nX * nY);

                // Compute x derivative second order central difference
                // f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / 12h

                int iPlus2  = grid.pbc_cell_coord<0>(i + 2);
                int iPlus1  = grid.pbc_cell_coord<0>(i + 1);
                int iMinus1 = grid.pbc_cell_coord<0>(i - 1);
                int iMinus2 = grid.pbc_cell_coord<0>(i - 2);

                real gradX = (         -UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, j, iPlus2) +
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, j, iPlus1) -
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, j, iMinus1) +
                                        UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, j, iMinus2));

                     gradX = gradX/(real(12.0)*h);

                int jPlus2  = grid.pbc_cell_coord<1>(j + 2);
                int jPlus1  = grid.pbc_cell_coord<1>(j + 1);
                int jMinus1 = grid.pbc_cell_coord<1>(j - 1);
                int jMinus2 = grid.pbc_cell_coord<1>(j - 2);

                real gradY = (         -UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, jPlus2, i) +
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, jPlus1, i) -
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, jMinus1, i) +
                                        UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, k, jMinus2, i));

                     gradY = gradY/(real(12.0)*h);

                int kPlus2  = grid.pbc_cell_coord<2>(k + 2);
                int kPlus1  = grid.pbc_cell_coord<2>(k + 1);
                int kMinus1 = grid.pbc_cell_coord<2>(k - 1);
                int kMinus2 = grid.pbc_cell_coord<2>(k - 2);

                real gradZ = (         -UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, kPlus2, j, i) +
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, kPlus1, j, i) -
                              real(8.0)*UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, kMinus1, j, i) +
                                        UAMMD_GET_3D_ROW_MAJOR(energy, 1, nY, nX, kMinus2, j, i));

                     gradZ = gradZ/(real(12.0)*h);

                // We store the -gradient (force)
                UAMMD_SET_3D_ROW_MAJOR(gradient, 1, nY, nX, k, j, i, -make_real3(gradX, gradY, gradZ));
            }
        }



    }

    struct HelixBoundaries_: public ExternalTabulated_{

        using StorageData = ExternalTabulated_::StorageData;

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data) {

            real helixPitch  = data.getParameter<real>("helixPitch");
            real helixRadius = data.getParameter<real>("helixRadius");
            real eps         = data.getParameter<real>("eps", 1.0);

            if(eps != 1.0 && eps != -1.0){
                System::log<System::CRITICAL>("[Helix Boundaries] "
                                              "eps must be either 1.0 or -1.0, but it is %f", eps);
            }

            System::log<System::MESSAGE>("[Helix Boundaries] "
                                         "HelixBoundaries potential with helixPitch = %f, helixRadius = %f, eps = %f",
                                         helixPitch, helixRadius, eps);

            real helixInnerRadius = data.getParameter<real>("helixInnerRadius");
            real K = data.getParameter<real>("K");

            System::log<System::MESSAGE>("[Helix Boundaries] "
                                         "Constraint potential with helixInnerRadius = %f, K = %f",
                                         helixInnerRadius, K);

            int nTurns = data.getParameter<int>("nTurns");
            real tMax = real(nTurns)*real(2*M_PI);

            System::log<System::MESSAGE>("[Helix Boundaries] "
                                         "Considering %d turns of the helix", nTurns);

            int nPointsHelix = data.getParameter<int>("nPointsHelix");

            System::log<System::MESSAGE>("[Helix Boundaries] "
                                         "Helix discretized with %d points", nPointsHelix);

            int nx = data.getParameter<int>("nx");
            int ny = data.getParameter<int>("ny");
            int nz = data.getParameter<int>("nz");

            Box box = gd->getEnsemble()->getBox();

            real dx = box.boxSize.x/real(nx);
            real dy = box.boxSize.y/real(ny);
            real dz = box.boxSize.z/real(nz);

            // Check all are equal
            if(dx != dy || dx != dz){
                System::log<System::CRITICAL>("[Helix Boundaries] The grid cells are not cubic."
                                              " dx = %f, dy = %f, dz = %f", dx, dy, dz);
            }

            real h = dx;

            // Check that the helix fits in the box
            if(helixPitch*nTurns != box.boxSize.z){
                System::log<System::CRITICAL>("[Helix Boundaries] "
                                              "Helix does not fit in the box."
                                              "helixPitch*nTurns = %f, box height = %f",
                                              helixPitch*nTurns, box.boxSize.z);
            }

            if(helixRadius > box.boxSize.x/real(2.0)){
                System::log<System::CRITICAL>("[Helix Boundaries] "
                                              "Helix radius is larger than half the box size."
                                              "helixRadius = %f, box width = %f",
                                              helixRadius, box.boxSize.x);
            }

            System::log<System::MESSAGE>("[Helix Boundaries] "
                                         "Starting grid boundary calculation");

            std::vector<real>  energy(nx*ny*nz, real(0.0));
            std::vector<real3> force(nx*ny*nz, make_real3(0.0, 0.0, 0.0));

            {
                thrust::device_vector<real>  energy_d(nx*ny*nz, real(0.0));
                thrust::device_vector<real3> force_d(nx*ny*nz, make_real3(0.0, 0.0, 0.0));


                real* energy_d_ptr = thrust::raw_pointer_cast(energy_d.data());

                int blockIterationsPerPoint = nPointsHelix/HelixBoundaries_ns::threadsPerBlock +
                                              int(nPointsHelix % HelixBoundaries_ns::threadsPerBlock != 0);

                dim3 gridDim(nx, ny, nz);
                HelixBoundaries_ns::fillGridDistanceEnergy2Helix<HelixBoundaries_ns::threadsPerBlock>
                <<<gridDim, HelixBoundaries_ns::threadsPerBlock>>>(helixPitch,
                                                                   helixRadius,
                                                                   eps,
                                                                   helixInnerRadius,
                                                                   K,
                                                                   tMax,
                                                                   box,
                                                                   h,
                                                                   blockIterationsPerPoint,
                                                                   nPointsHelix,
                                                                   energy_d_ptr);

                cudaDeviceSynchronize();

                Grid grid(box, h);

                int Nthreads = 512;
                int Nblocks  = (nx*ny*nz)/Nthreads + (((nx*ny*nz)%Nthreads)?1:0);

                real3* force_d_ptr = thrust::raw_pointer_cast(force_d.data());

                HelixBoundaries_ns::computeForce<<<Nblocks, Nthreads>>>(box,
                                                                        grid,
                                                                        energy_d_ptr,
                                                                        force_d_ptr);

                cudaDeviceSynchronize();
                thrust::copy(energy_d.begin(), energy_d.end(), energy.begin());
                thrust::copy(force_d.begin(), force_d.end(), force.begin());
            }

            std::vector<int> id_x(nx*ny*nz, 0);
            std::vector<int> id_y(nx*ny*nz, 0);
            std::vector<int> id_z(nx*ny*nz, 0);

            std::vector<real>  energy_rearranged(nx*ny*nz, real(0.0));
            std::vector<real3> force_rearranged(nx*ny*nz, make_real3(0.0, 0.0, 0.0));

            int index=0;
            for(int k = 0; k < nz; k++){
                for(int j = 0; j < ny; j++){
                    for(int i = 0; i < nx; i++){
                        id_x[index] = i;
                        id_y[index] = j;
                        id_z[index] = k;
                        energy_rearranged[index] = UAMMD_GET_3D_ROW_MAJOR(energy.data(), nz, ny, nx, k, j, i);
                        force_rearranged[index]  = UAMMD_GET_3D_ROW_MAJOR(force.data(), nz, ny, nx, k, j, i);
                        index++;
                    }
                }
            }

            return ExternalTabulated_::getStorageData(gd, pg,
                                                      nx, ny, nz,
                                                      1.0,
                                                      id_x,id_y,id_z,
                                                      energy_rearranged, force_rearranged);
        };

        struct ComputationalData : public ExternalTabulated_::ComputationalData {
            real lambda;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            static_cast<ExternalTabulated_::ComputationalData&>(computational) =
            ExternalTabulated_::getComputationalData(gd, pg, storage, computables, st);

            // Get lambda
            computational.lambda = gd->getEnsemble()->getLambda();

            return computational;
        };

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
                real e = ExternalTabulated_::energy(index_i, computational);
                     e = e*computational.lambda*computational.lambda;
                return e;
        }

        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
                real3 f = ExternalTabulated_::force(index_i, computational);
                      f = f*computational.lambda*computational.lambda;
                return f;
        }

        static inline __device__ real lambdaDerivative(int index_i,const ComputationalData& computational){
                real e  = ExternalTabulated_::energy(index_i, computational);
                real ld = real(2.0)*computational.lambda*e;
                return ld;
        }

    };

    using HelixBoundaries = ExternalLambda_<HelixBoundaries_>;

}}}}

REGISTER_SINGLE_INTERACTOR(
    External,HelixBoundaries,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::HelixBoundaries>
)
