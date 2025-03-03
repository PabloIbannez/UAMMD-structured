#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"

#include "Interactor/Interactor.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Matrices.cuh"

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
            thrust::device_vector<real4> buffer_d;

            int pixelBufferSize;
            thrust::device_vector<real2> pixelPos_d;
            thrust::device_vector<real>  pixelValue_d;

            //

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

                //

            }

            ~AFMimageInteractor(){}

            void sum(Computables comp,cudaStream_t st) override {

                if(comp.energy == true || comp.force){
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
