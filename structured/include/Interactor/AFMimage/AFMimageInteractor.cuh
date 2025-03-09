#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"

#include "Interactor/Interactor.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Matrices.cuh"

#include "Utils/Maths/MatrixIndexing.cuh"

#include "cub/device/device_segmented_reduce.cuh"

//TODO: Ensure we work with ids not indices !!!! also for the buffer !!!! <--- IMPORTANT (Check populateBuffer function !!!!)
//TODO: I think the Npx and Npy should be extended to match the box size

namespace uammd{
namespace structured{
namespace Interactor{

    namespace AFMimageInteractor_ns {

        __device__ real phi(const real3& partPos, const real2& pixelPos,
                            const real& sigma, const real& gamma, const real& tipRadius,
                            Box box){

            const real2 rpartpix  = make_real2(box.apply_pbc(make_real3(partPos.x,partPos.y,0.0)-make_real3(pixelPos,0.0)));
            const real2 rpartpix2 = rpartpix*rpartpix;

            const real xy = exp(-0.5*(rpartpix2.x + rpartpix2.y)/(sigma*sigma));
            const real z  = exp((partPos.z + tipRadius)/gamma);

            return xy*z;
        }

        __device__ real3 phi_derivative(const real3& partPos, const real2& pixelPos,
                                        const real& sigma, const real& gamma, const real& tipRadius,
                                        Box box){
            //Compute the derivative of phi with respect to the particle position
            const real2 rpartpix  = make_real2(box.apply_pbc(make_real3(partPos.x,partPos.y,0.0)-make_real3(pixelPos,0.0)));
            const real2 rpartpix2 = rpartpix*rpartpix;

            const real xy = exp(-0.5*(rpartpix2.x + rpartpix2.y)/(sigma*sigma));
            const real z  = exp((partPos.z + tipRadius)/gamma);

            const real phi_val = xy*z;

            const real dx = -(rpartpix.x/(sigma*sigma))*phi_val;
            const real dy = -(rpartpix.y/(sigma*sigma))*phi_val;
            const real dz = (1.0/gamma)*phi_val;

            return make_real3(dx,dy,dz);
        }

        __device__ real image(const real& gamma, const real& pixelSum){
            return gamma*log(real(1.0) + pixelSum);
        }

        __device__ real3 image_derivative(const real& gamma, const real3& phi_derivative_val, const real& pixelSum){
            return gamma*phi_derivative_val/(real(1.0) + pixelSum);
        }

        __global__ void populateBuffer(real4* buffer,
                                       real4* pos, int* batchId,
                                       real sigma, real gamma, real tipRadius,
                                       real2* pixelPos, real* pixelValue,
                                       int Npx, int Npy, int N, int Nbatches, Box box){
            // This kernel is executed by Npx*Npy*N threads

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int n  = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && n < N) {

                real3 partPos = make_real3(pos[n]);
                int   currentBatch = batchId[n];

                real2 pPos = UAMMD_GET_3D_ROW_MAJOR(pixelPos,Nbatches,Npx,Npy,currentBatch,px,py);

                real  phi_val = phi(partPos,pPos,sigma,gamma,tipRadius,box);
                real3 phi_derivative_val = phi_derivative(partPos,pPos,sigma,gamma,tipRadius,box);

                real4 phi_buffer = make_real4(phi_derivative_val,phi_val);

                UAMMD_SET_3D_ROW_MAJOR(buffer,N,Npx,Npy,n,px,py,phi_buffer);
            }

        }

        __global__ void updateBuffer(real4* buffer,
                                     int*   batchId,
                                     real*  pixelValueTheoretical, // It contains sum_(i=1)^N_i phi_i for each pixel
                                     real   gamma,
                                     int Npx, int Npy, int N, int Nbatches){
            // This kernel is executed by Npx*Npy*N threads

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int n  = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && n < N) {

                int currentBatch = batchId[n];

                real4 buffer_val = UAMMD_GET_3D_ROW_MAJOR(buffer,N,Npx,Npy,n,px,py);
                real pixelSum    = UAMMD_GET_3D_ROW_MAJOR(pixelValueTheoretical,Nbatches,Npx,Npy,currentBatch,px,py);

                real3 phi_derivative_val = make_real3(buffer_val);
                real3 image_derivative_val = image_derivative(gamma,phi_derivative_val,pixelSum);

                // Update the buffer
                real4 updated_buffer_val = make_real4(image_derivative_val,buffer_val.w);

                UAMMD_SET_3D_ROW_MAJOR(buffer,N,Npx,Npy,n,px,py,updated_buffer_val);
            }

        }

        // toReduceValue = buffer[reductionIds[i]].w
        __global__ void prepareToReduceValueKernel(real4* buffer,
                                                   int* reductionIds,
                                                   real* toReduceValue,
                                                   int N){
            const int i = blockIdx.x * blockDim.x + threadIdx.x;

            if(i < N){
                int index = reductionIds[i];
                toReduceValue[i] = buffer[index].w;
            }
        }

        __global__ void computeImage(real* pixelValueTheoretical,
                                     real gamma,
                                     int Npx, int Npy, int Nbatches){

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int nb = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && nb < Nbatches) { //Note we use Nbatches here, instead of N
                real pixelSum = UAMMD_GET_3D_ROW_MAJOR(pixelValueTheoretical,Nbatches,Npx,Npy,nb,px,py);
                real pixelValue = image(gamma,pixelSum);
                UAMMD_SET_3D_ROW_MAJOR(pixelValueTheoretical,Nbatches,Npx,Npy,nb,px,py,pixelValue); // Note we read and write from the same array. No problem here.
            }
        }

    }

    template<class AFMimageType>
    class AFMimageInteractor: public Interactor{

        private:

            const int THREADS_PER_X = 8;
            const int THREADS_PER_Y = 8;
            const int THREADS_PER_N = 8;

            int N;

            int Npx;
            int Npy;
            int Nbatches;

            real sigma;
            real gamma;
            real tipRadius;

            std::map<int,real2> batchOffset;

            std::map<int,real2> batchMin;
            std::map<int,real2> batchMax;

            std::map<int,real2> batchResolution;

            std::map<int,std::map<std::pair<int,int>,real3>> batchId2pixels;

            std::shared_ptr<GlobalData>   gd;
            std::shared_ptr<AFMimageType> afm;

            // GPU data

            int bufferSize;
            // buffer stores the values of the phi function and its derivative
            // The size of the buffer is (N,Npx,Npy) and it stores for each particle the values of the phi function and its derivative
            // Note that here it uses the global index of the particle, meaning that the buffer ignores the batch information itself
            thrust::device_vector<real4> buffer_d;

            int pixelBufferSize;
            // The two following vectors are constant.
            // They stored the input data in a 3D array (batchId,x,y) to be used in the kernel
            //
            // The first vector encoded a 3d matrix of real2 (x,y).
            // The size of the matrix is (Nbatches,Npx,Npy) and it stores for each batch the x and y values of the pixels
            // Meaning, it transform from integer coordinates, i,j, to real2 coordinates, x,y for each batch
            //
            // The second vector encodes a 3d matrix of real (height).
            // The size of the matrix is (Nbatches,Npx,Npy) and it stores for each batch the height of the pixels
            // Meaning, it transform from integer coordinates, i,j, to real height for each batch
            thrust::device_vector<real2> pixelPos_d;
            thrust::device_vector<real>  pixelValue_d;

            thrust::device_vector<real>  pixelValueTheoretical_d;

            thrust::device_vector<int> reductionIds_d;
            thrust::device_vector<int> reductionIdsSegmentStart_d;

            // cub segmented reduction stuff
            // int cub_num_segments = Npx*Npy*Nbatches;
            size_t cub_temp_storage_bytes = 0;
            void* cub_temp_storage = nullptr;

            thrust::device_vector<real> toReduceValue_d;

            void prepareToReduceValue(cudaStream_t st){

                int Nthreads = 256;
                int Nblocks  = bufferSize/Nthreads + ((bufferSize%Nthreads)?1:0);

                real4* buffer_ptr        = thrust::raw_pointer_cast(buffer_d.data());
                int*   reductionIds_ptr  = thrust::raw_pointer_cast(reductionIds_d.data());
                real*  toReduceValue_ptr = thrust::raw_pointer_cast(toReduceValue_d.data());

                AFMimageInteractor_ns::prepareToReduceValueKernel<<<Nblocks,Nthreads,0,st>>>(buffer_ptr,
                                                                                             reductionIds_ptr,
                                                                                             toReduceValue_ptr,
                                                                                             bufferSize);
            }

            /////////////////////////////

            bool warningEnergy = false;
            bool warningForce  = false;

        public:

            AFMimageInteractor(std::shared_ptr<GlobalData>    gd,
                               std::shared_ptr<ParticleGroup> pg,
                               DataEntry& data,
                               std::string name):Interactor(pg,"AFMimageInteractor: \"" +name+"\""),
                                                 gd(gd){

                N = pg->getParticleData()->getNumParticles();
                if(N != pg->getNumberParticles()){
                    System::log<System::CRITICAL>("[AFMimageInteractor] (%s) The current implementation of AFMimageInteractor only supports"
                                                  " ParticleGroup that contains all the particles in the simulation",
                                                  name.c_str());
                }

                Npx = data.getParameter<int>("Npx");
                Npy = data.getParameter<int>("Npy");

                sigma = data.getParameter<real>("sigma");
                gamma = data.getParameter<real>("gamma");
                tipRadius = data.getParameter<real>("tipRadius");

                System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Initialized with Npx=%d, Npy=%d",
                                             name.c_str(),Npx,Npy);

                std::vector<real> px      = data.getData<real>("x");
                std::vector<real> py      = data.getData<real>("y");
                std::vector<real> height  = data.getData<real>("height");

                // Start batchId checking
                std::vector<int> batchIds;
                if(data.isDataAdded("batchId")){

                    batchIds = data.getData<int>("batchId");

                    std::vector<int> uniqueBatchIds = batchIds;
                    std::sort(uniqueBatchIds.begin(),uniqueBatchIds.end());
                    uniqueBatchIds.erase(std::unique(uniqueBatchIds.begin(),uniqueBatchIds.end()),
                                         uniqueBatchIds.end());

                    // Check that batchIds start from 0 and are consecutive
                    if(uniqueBatchIds[0] != 0){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) BatchIds must start from 0",
                                                      name.c_str());
                    }
                    for(int i=1;i<uniqueBatchIds.size();i++){
                        if(uniqueBatchIds[i] != uniqueBatchIds[i-1]+1){
                            System::log<System::CRITICAL>("[AFMimageInteractor] (%s) BatchIds must be consecutive",
                                                          name.c_str());
                        }
                    }
                    Nbatches = uniqueBatchIds.size();

                    // Check that the number of batches is consistent with the number of batchIds
                    auto batId = pd->getBatchId(access::location::cpu,access::mode::read);
                    std::vector<int> uniqueBatchIds_pd(batId.begin(),batId.end());
                    std::sort(uniqueBatchIds_pd.begin(),uniqueBatchIds_pd.end());
                    uniqueBatchIds_pd.erase(std::unique(uniqueBatchIds_pd.begin(),uniqueBatchIds_pd.end()),
                                         uniqueBatchIds_pd.end());

                    // uniqueBatchIds and uniqueBatchIds_pd must be the same
                    if(uniqueBatchIds.size() != uniqueBatchIds_pd.size()){
                        int nBatchesData = uniqueBatchIds.size();
                        int nBatchesPD   = uniqueBatchIds_pd.size();
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Number of batches is inconsistent with batchIds (Data: %d, PD: %d)",
                                                      name.c_str(),nBatchesData,nBatchesPD);
                        if(uniqueBatchIds != uniqueBatchIds_pd){
                            System::log<System::CRITICAL>("[AFMimageInteractor] (%s) BatchIds are inconsistent",
                                                          name.c_str());
                        }
                    }

                } else {

                    Nbatches = 1;
                    System::log<System::MESSAGE>("[AFMimageInteractor] (%s) No batch information found",
                                                 name.c_str());

                    auto batId = pd->getBatchId(access::location::cpu,access::mode::read);
                    for(int i=0;i<batId.size();i++){
                        // Check if all batchIds are 0
                        if(batId[i] != 0){
                            System::log<System::CRITICAL>("[AFMimageInteractor] (%s) All batchIds must be 0",
                                                          name.c_str());
                        }
                    }

                    batchIds.resize(px.size(),0);
                }

                System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Found %d batches",
                                             name.c_str(),Nbatches);

                // Since we are assuming that each batch has the same number of pixels, Npx*Npy,
                // we can check that the number of pixels is consistent with the number of batches

                if(px.size() != Npx*Npy*Nbatches || py.size() != Npx*Npy*Nbatches || height.size() != Npx*Npy*Nbatches){
                    System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Number of pixels is inconsistent with the given image size (Npx=%d, Npy=%d, total=%d) "
                                                  "and the number of batches (Nbatches=%d). The total given data should be Npx*Npy*Nbatches=%d",
                                                  name.c_str(),Npx,Npy,Npx*Npy,Nbatches,Npx*Npy*Nbatches);
                }

                // At this point batchIds is consistent with the batch information in the ParticleData
                // and it is loaded in the vector batchIds


                // Check pixels
                // We create a data structure that associates each batch with another map that associates
                // each pixel (using int coordinates) with the pixel value (real3, x,y,height)
                // std::map<int,std::map<std::pair<int,int>,real3>> batchId2pixels;

                // Initialize the data structure
                for(int i=0;i<Nbatches;i++){
                    batchId2pixels[i] = std::map<std::pair<int,int>,real3>();
                }

                for(int currentBatch=0;currentBatch<Nbatches;currentBatch++){
                    std::vector<real3> batchPixels(Npx*Npy);
                    std::set<real> xUniqueValues;
                    std::set<real> yUniqueValues;
                    for(int i=0;i<Npx*Npy*Nbatches;i++){
                        if(batchIds[i] == currentBatch){
                            batchPixels.push_back(make_real3(px[i],py[i],height[i]));
                            xUniqueValues.insert(px[i]);
                            yUniqueValues.insert(py[i]);
                        }
                    }
                    // Check that the number of unique x and y values is consistent with Npx and Npy
                    if(xUniqueValues.size() != Npx || yUniqueValues.size() != Npy){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Number of given x and y values is inconsistent with the given image size (Npx=%d, Npy=%d)",
                                                      name.c_str(),Npx,Npy);
                    }

                    // Sort the unique values
                    std::vector<real> xUniqueValuesSorted(xUniqueValues.begin(),xUniqueValues.end());
                    std::vector<real> yUniqueValuesSorted(yUniqueValues.begin(),yUniqueValues.end());
                    std::sort(xUniqueValuesSorted.begin(),xUniqueValuesSorted.end());
                    std::sort(yUniqueValuesSorted.begin(),yUniqueValuesSorted.end());

                    // Determine the offset and resolution of the batch
                    real2 offset = make_real2(xUniqueValuesSorted[0],yUniqueValuesSorted[0]);
                    batchOffset[currentBatch] = offset;

                    real2 resolution = make_real2(xUniqueValuesSorted[1]-xUniqueValuesSorted[0],
                                                  yUniqueValuesSorted[1]-yUniqueValuesSorted[0]);
                    // Check that the resolution is consistent across all the pixels pairs
                    for(int i=1;i<Npx;i++){
                        if(xUniqueValuesSorted[i]-xUniqueValuesSorted[i-1] != resolution.x){
                            System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Inconsistent resolution in x for batch %d",
                                                          name.c_str(),currentBatch);
                        }
                    }
                    for(int i=1;i<Npy;i++){
                        if(yUniqueValuesSorted[i]-yUniqueValuesSorted[i-1] != resolution.y){
                            System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Inconsistent resolution in y for batch %d",
                                                          name.c_str(),currentBatch);
                        }
                    }
                    // Check image is inside the box
                    Box box = gd->getEnsemble()->getBox();
                    // We assume that the position of the pixel is the center of the pixel
                    real2 pixelSize = make_real2(resolution.x/2.0,resolution.y/2.0);
                    real2 pixelMin  = offset - pixelSize;
                    real2 pixelMax  = offset + pixelSize + make_real2(resolution.x*(Npx-1),resolution.y*(Npy-1));

                    real boxMinX = -box.boxSize.x/2.0;
                    real boxMaxX =  box.boxSize.x/2.0;
                    real boxMinY = -box.boxSize.y/2.0;
                    real boxMaxY =  box.boxSize.y/2.0;

                    if(pixelMin.x < boxMinX || pixelMax.x > boxMaxX || pixelMin.y < boxMinY || pixelMax.y > boxMaxY){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Image is outside the box for batch %d",
                                                      name.c_str(),currentBatch);
                    }

                    batchMin[currentBatch] = pixelMin;
                    batchMax[currentBatch] = pixelMax;

                    // Everything is consistent, we can store the resolution
                    batchResolution[currentBatch] = resolution;

                    System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Batch %d has resolution (%f,%f)",
                                                 name.c_str(),currentBatch,resolution.x,resolution.y);

                    // Populate the data structure
                    for(int i=0;i<Npx*Npy;i++){
                        int x = (px[i]-offset.x)/resolution.x;
                        int y = (py[i]-offset.y)/resolution.y;
                        std::pair<int,int> coord = std::make_pair(x,y);
                        batchId2pixels[currentBatch][coord] = batchPixels[i];
                    }

                }

                // Prepare GPU data
                bufferSize = Npx*Npy*N;
                buffer_d.resize(bufferSize);

                pixelBufferSize = Npx*Npy*Nbatches;
                thrust::host_vector<real2> pixelPos_h(pixelBufferSize);
                thrust::host_vector<real>  pixelValue_h(pixelBufferSize);

                // Both pixelPos_d and pixelValue_d are constant, we can populate them here
                for(int currentBatch=0;currentBatch<Nbatches;currentBatch++){
                    real2 offset     = batchOffset[currentBatch];
                    real2 resolution = batchResolution[currentBatch];
                    for(int i=0;i<Npx;i++){
                    for(int j=0;j<Npy;j++){
                        // #define UAMMD_SET_3D_ROW_MAJOR(matrix, depth, rows, cols, dep, row, col, value) (matrix)[(dep) * (rows) * (cols) + (row) * (cols) + (col)] = (value)
                        real3 pixel = batchId2pixels[currentBatch][std::make_pair(i,j)];

                        real2 pixelPos;
                        pixelPos.x = pixel.x;
                        pixelPos.y = pixel.y;
                        UAMMD_SET_3D_ROW_MAJOR(pixelPos_h,Nbatches,Npx,Npy,currentBatch,i,j,pixelPos);

                        real  pixelValue = pixel.z;
                        UAMMD_SET_3D_ROW_MAJOR(pixelValue_h,Nbatches,Npx,Npy,currentBatch,i,j,pixelValue);
                    }}
                }

                pixelPos_d   = pixelPos_h;
                pixelValue_d = pixelValue_h;

                // Let's prepare the indexing arrays for reductions
                // Mainly we have to reduce the columns of the buffer
                // and store the ruduced values in the array pixelValueTheoretical_d
                // pixelValueTheoretical_d has size (Nbatches,Npx,Npy) (as pixelValue_d)
                // SO we have to determine the indices of the buffer that correspond to each pixel.
                // For the pixel (currentBatch,px,py) the indices in the buffer that their values
                // reduction give the pixel value.

                // First, we initialize the result array. The array for the theoretical image
                pixelValueTheoretical_d.resize(Nbatches*Npx*Npy);

                // Now we are going to create an auxiliary array that associates each batch with the
                // ids of the particles that belong to that batch.

                std::map<int,std::vector<int>> batchId2ids;
                for(int i=0;i<Nbatches;i++){
                    batchId2ids[i] = std::vector<int>();
                }

                auto ids     = this->pd->getId(access::location::cpu, access::mode::read);
                auto batchId = this->pd->getBatchId(access::location::cpu, access::mode::read);

                for(int i=0;i<N;i++){
                    batchId2ids[batchId[i]].push_back(ids[i]);
                }

                // The array reductionIds_h will store the indices of the buffer that correspond to each pixel
                // The array reductionIdsSegmentStart_h will store the starting index of each segment in the reductionIds_h
                thrust::host_vector<int> reductionIds_h(N*Npx*Npy,-1);
                thrust::host_vector<int> reductionIdsSegmentStart_h(Nbatches*Npx*Npy,-1);

                int previousId  = 0;
                int currentSegmentStart = 0;

                // First we iterate (IN ORDER) over the elements of pixelValueTheoretical_d
                int pixelValueTheoreticalSize = Npx*Npy*Nbatches;
                for(int pixelValueTheoretical_index=0;
                        pixelValueTheoretical_index<pixelValueTheoreticalSize;
                        pixelValueTheoretical_index++){

                    int3 coordinates = MatrixIndexing::rowMajor3Dcoordinates(Nbatches,Npx,Npy,pixelValueTheoretical_index);

                    int currentBatch = coordinates.x;
                    int px = coordinates.y;
                    int py = coordinates.z;

                    // Now we have to determine the indices in the buffer that correspond to the pixel (currentBatch,px,py)
                    for(int i=0;i<batchId2ids[currentBatch].size();i++){
                        int id = batchId2ids[currentBatch][i];
                        int index = MatrixIndexing::rowMajor3Dindex(N,Npx,Npy,id,px,py);
                        reductionIds_h[previousId] = index;
                        previousId++;
                    }
                    reductionIdsSegmentStart_h[pixelValueTheoretical_index] = currentSegmentStart;
                    currentSegmentStart += batchId2ids[currentBatch].size();
                }

                // Check if all elements of reductionIds_h and reductionIdsSegmentStart_h have been filled
                for(int i=0;i<reductionIds_h.size();i++){
                    if(reductionIds_h[i] == -1){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Reduction indices have not been filled",
                                                      name.c_str());
                    }
                }
                for(int i=0;i<reductionIdsSegmentStart_h.size();i++){
                    if(reductionIdsSegmentStart_h[i] == -1){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Reduction indices segment start have not been filled",
                                                      name.c_str());
                    }
                }

                reductionIds_d    = reductionIds_h;
                reductionIdsSegmentStart_d = reductionIdsSegmentStart_h;

                System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Reduction indices prepared",
                                             name.c_str());
                // First we reserve memory for the auxiliary array to be used in the reduction
                toReduceValue_d.resize(N*Npx*Npy);

                // At this point we have to prepare all the stuff needed by cub for the reduction
                // We use DeviceSegmentedReduce for the reduction
                {
                    prepareToReduceValue(0); // This is not necessary here, but it is a good way to check if the kernel works

                    real* toReduceValue_ptr = thrust::raw_pointer_cast(toReduceValue_d.data());
                    real* pixelValueTheoretical_ptr = thrust::raw_pointer_cast(pixelValueTheoretical_d.data());

                    // Now we have to determine the size of the temporary storage needed by cub
                    cub::DeviceSegmentedReduce::Sum(cub_temp_storage,cub_temp_storage_bytes,
                                                    toReduceValue_ptr,pixelValueTheoretical_ptr, //d_in, d_out
                                                    Npx*Npy*Nbatches, // num_segments
                                                    reductionIdsSegmentStart_d.begin(),
                                                    reductionIdsSegmentStart_d.begin()+1,
                                                    0); // stream
                    cudaDeviceSynchronize();

                    cudaMalloc(&cub_temp_storage,cub_temp_storage_bytes*sizeof(real));

                }

            }

            ~AFMimageInteractor(){
                if(cub_temp_storage != nullptr){
                    cudaFree(cub_temp_storage);
                }
            }

            void sum(Computables comp,cudaStream_t st) override {

                if(comp.energy == true || comp.force == true){
                    // Common for both energy and force
                    System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Computing common data",
                                                 name.c_str());

                    // Populate the buffer
                    dim3 threadsPerBlock(THREADS_PER_X,THREADS_PER_Y,THREADS_PER_N);
                    dim3 numBlocks;

                    numBlocks.x = Npx/threadsPerBlock.x + ((Npx%threadsPerBlock.x)?1:0);
                    numBlocks.y = Npy/threadsPerBlock.y + ((Npy%threadsPerBlock.y)?1:0);
                    numBlocks.z = N/threadsPerBlock.z   + ((N%threadsPerBlock.z)?1:0);

                    real4* pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite).raw();
                    int* batchId = this->pd->getBatchId(access::location::gpu, access::mode::readwrite).raw();

                    real4* buffer = thrust::raw_pointer_cast(buffer_d.data());

                    real2* pixelPos = thrust::raw_pointer_cast(pixelPos_d.data());
                    real*  pixelValue = thrust::raw_pointer_cast(pixelValue_d.data());

                    AFMimageInteractor_ns::populateBuffer<<<numBlocks,threadsPerBlock,0,st>>>(buffer,
                                                                                              pos,batchId,
                                                                                              sigma,gamma,tipRadius,
                                                                                              pixelPos,pixelValue,
                                                                                              Npx,Npy,N,Nbatches,
                                                                                              gd->getEnsemble()->getBox());
                    {
                        prepareToReduceValue(st);

                        real* toReduceValue_ptr = thrust::raw_pointer_cast(toReduceValue_d.data());
                        real* pixelValueTheoretical_ptr = thrust::raw_pointer_cast(pixelValueTheoretical_d.data());

                        // Now we have to determine the size of the temporary storage needed by cub
                        cub::DeviceSegmentedReduce::Sum(cub_temp_storage,cub_temp_storage_bytes,
                                                        toReduceValue_ptr,pixelValueTheoretical_ptr, //d_in, d_out
                                                        Npx*Npy*Nbatches, // num_segments
                                                        reductionIdsSegmentStart_d.begin(),
                                                        reductionIdsSegmentStart_d.begin()+1,
                                                        st); // stream
                        cudaDeviceSynchronize();
                    }

                    // At this point we have the buffer with the values of the phi function and its derivative
                    // and we have the quantity sum_(i=1)^N_i phi_i for each pixel stored in pixelValueTheoretical_d

                    // Now we have to update the buffer (we use the same configuration as before)
                    AFMimageInteractor_ns::updateBuffer<<<numBlocks,threadsPerBlock,0,st>>>(buffer,
                                                                                            batchId,
                                                                                            thrust::raw_pointer_cast(pixelValueTheoretical_d.data()),
                                                                                            gamma,
                                                                                            Npx,Npy,N,Nbatches);

                    // Now each element of the buffer has the value of phi, and the derivative of phi respect to the particle position

                    // Now we have to compute the image
                    numBlocks.z = Nbatches/threadsPerBlock.z + ((Nbatches%threadsPerBlock.z)?1:0);

                    real* pixelValueTheoretical_ptr = thrust::raw_pointer_cast(pixelValueTheoretical_d.data());

                    AFMimageInteractor_ns::computeImage<<<numBlocks,threadsPerBlock,0,st>>>(pixelValueTheoretical_ptr,
                                                                                            gamma,
                                                                                            Npx,Npy,Nbatches);
                    // At this point we have the theoretical image stored in pixelValueTheoretical_d

                    // Now we have to compute the comparison and the normalization operations. These operations can vary
                    // depending on the type of the potential.
                }

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<AFMimageType>::value){

                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[AFMimageInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }
                }

                if(comp.force == true){

                    if constexpr (has_getForceTransverser<AFMimageType>::value){

                    } else {
                        if(!warningForce){
                            System::log<System::WARNING>("[AFMimageInteractor] (%s) Requested non-implemented transverser (force)",
                                                         name.c_str());
                            warningForce = true;
                        }
                    }
                }
            }
    };

}}}
