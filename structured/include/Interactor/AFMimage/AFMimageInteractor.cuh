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

//TODO: It should be possible to perform the reduction of pixelImage using only the component w. But I have not been able to explain this to cub
//TODO: I think the Npx and Npy should be extended to match the box size
//TODO: Divide the energy by Nparticles to get the average energy per particle. This has to be done in AFMimage.cuh

namespace uammd{
namespace structured{
namespace Interactor{

    namespace AFMimageInteractor_ns {

        struct Real4Sum {
            __host__ __device__
                real4 operator()(const real4 &a, const real4 &b) const {
                    return a + b;
                }
        };

        template<class ComputationalData, class AFMimageTransverser>
        __global__ void transverseAFMimage(const ComputationalData computational,
                                           AFMimageTransverser afmImageTransverser,
                                           int* ids,
                                           real4* particleImageReductionResult,
                                           real4* particleImageReductionResultAux,
                                           int N) {
            // This kernel is executed by N threads
            const int index = blockIdx.x * blockDim.x + threadIdx.x;

            if (index < N) {
                const int id = ids[index];

                const real4 reductionResult    = particleImageReductionResult[id];
                const real4 reductionResultAux = particleImageReductionResultAux[id];

                const real opGen         = reductionResult.w;
                const real opRefGen      = reductionResultAux.w;
                const real3 opRefGenDer  = make_real3(reductionResultAux);
                const real3 opGenGenDer  = make_real3(reductionResult);

                typename AFMimageTransverser::resultType quantity = afmImageTransverser.zero();
                afmImageTransverser.accumulate(quantity,afmImageTransverser.compute(index,computational,
                                                                                    opGen,opRefGen,
                                                                                    opRefGenDer,opGenGenDer));

                afmImageTransverser.set(index,quantity);
            }

        }

        __device__ real phi(const real3& partPos, const real2& pixelPos,
                            const real& sigma, const real& gamma, const real& tipRadius,
                            Box box){

            const real2 rpartpix  = make_real2(box.apply_pbc(make_real3(partPos.x,partPos.y,0.0)-make_real3(pixelPos,0.0)));
            const real2 rpartpix2 = rpartpix*rpartpix;

            const real xy = exp(-real(0.5)*(rpartpix2.x + rpartpix2.y)/(sigma*sigma));
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
                                       int* ids,
                                       real4* pos, int* batchId,
                                       real sigma, real gamma, real tipRadius,
                                       real2* pixelPos,
                                       int Npx, int Npy, int N, int Nbatches, Box box){
            // This kernel is executed by Npx*Npy*N threads

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int n  = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && n < N) {

                int   currentId = ids[n];

                real3 partPos = make_real3(pos[n]);
                int   currentBatch = batchId[n];

                real2 pPos = UAMMD_GET_3D_ROW_MAJOR(pixelPos,Nbatches,Npx,Npy,currentBatch,px,py);

                real  phi_val = phi(partPos,pPos,sigma,gamma,tipRadius,box);
                real3 phi_derivative_val = phi_derivative(partPos,pPos,sigma,gamma,tipRadius,box);

                real4 phi_buffer = make_real4(phi_derivative_val,phi_val);

                UAMMD_SET_3D_ROW_MAJOR(buffer,N,Npx,Npy,currentId,px,py,phi_buffer);
            }

        }

        template<class AFMimageType>
        __global__ void updateBuffer(real4* buffer,
                                     real4* bufferAux,
                                     int* ids,
                                     real4* pixelReductionResult, // It contains sum_(i=1)^N_i phi_i for each pixel at .w
                                     real*  referenceImages, // It contains the reference images for each batch
                                     int* batchId,
                                     real   gamma,
                                     int Npx, int Npy, int N, int Nbatches){
            // This kernel is executed by Npx*Npy*N threads

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int n  = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && n < N) {

                int   currentId = ids[n];
                int currentBatch = batchId[n];

                real4 buffer_val = UAMMD_GET_3D_ROW_MAJOR(buffer,N,Npx,Npy,currentId,px,py);
                real4 pixelSum_ = UAMMD_GET_3D_ROW_MAJOR(pixelReductionResult,Nbatches,Npx,Npy,currentBatch,px,py);
                real  pixelSum = pixelSum_.w;

                real  gen_image_val     = image(gamma,pixelSum);

                real3 phi_derivative_val   = make_real3(buffer_val);
                real3 gen_image_derivative_val = image_derivative(gamma,phi_derivative_val,pixelSum);

                real  ref_image_val = UAMMD_GET_3D_ROW_MAJOR(referenceImages,Nbatches,Npx,Npy,currentBatch,px,py);

                // Update the buffer
                real GEN     = AFMimageType::OperationOverGenerated(gen_image_val);
                real REF_GEN = AFMimageType::OperationOverReferenceAndGenerated(ref_image_val,gen_image_val);
                real3 REF_GEN_derivative = AFMimageType::OperationOverReferenceAndGeneratedDerivative(ref_image_val,gen_image_derivative_val);
                real3 GEN_GEN_derivative = AFMimageType::OperationOverGeneratedAndGeneratedDerivative(gen_image_val,gen_image_derivative_val);

                real4 result1 = make_real4(GEN_GEN_derivative,GEN);
                real4 result2 = make_real4(REF_GEN_derivative,REF_GEN);

                UAMMD_SET_3D_ROW_MAJOR(buffer,N,Npx,Npy,currentId,px,py,result1);
                UAMMD_SET_3D_ROW_MAJOR(bufferAux,N,Npx,Npy,currentId,px,py,result2);

            }

        }

        __global__ void computeGeneratedImagesFromPixelReductionResults(real4* pixelReductionResult,
                                                                        real* generatedImages,
                                                                        real gamma,
                                                                        int Npx, int Npy, int Nbatches){

            // Compute the indices for each dimension.
            const int px = blockIdx.x * blockDim.x + threadIdx.x;
            const int py = blockIdx.y * blockDim.y + threadIdx.y;
            const int nb = blockIdx.z * blockDim.z + threadIdx.z;

            // Check if indices are within the matrix dimensions.
            if (px < Npx && py < Npy && nb < Nbatches) { //Note we use Nbatches here, instead of N
                real4 pixelSum = UAMMD_GET_3D_ROW_MAJOR(pixelReductionResult,Nbatches,Npx,Npy,nb,px,py);
                real pixelValue = image(gamma,pixelSum.w); // Note we use component w of pixelSum !!!. This is the sum of all the phi values
                UAMMD_SET_3D_ROW_MAJOR(generatedImages,Nbatches,Npx,Npy,nb,px,py,pixelValue);
            }
        }

        class AFMimage{

            private:

                int Npx,Npy,Nbatches;

                std::vector<real3> xy_height;

            public:

                AFMimage(const thrust::host_vector<real2>& pixelPos,
                         const thrust::host_vector<real>&  image,
                         int Npx, int Npy, int Nbatches):Npx(Npx),Npy(Npy),Nbatches(Nbatches){

                    xy_height.resize(Npx*Npy*Nbatches);

                    const real2* pixelPos_ptr = thrust::raw_pointer_cast(pixelPos.data());
                    const real*  image_ptr    = thrust::raw_pointer_cast(image.data());

                    for(int i=0;i<Npx*Npy*Nbatches;i++){
                        real2 ppos = pixelPos[i];
                        xy_height[i] = make_real3(ppos.x,ppos.y,image[i]);
                    }

                }

                real3 operator()(int px, int py, int batch){
                    if(px < 0 || px >= Npx || py < 0 || py >= Npy || batch < 0 || batch >= Nbatches){
                        System::log<System::CRITICAL>("[AFMimageInteractor::image] Error in the indices (px=%d, py=%d, batch=%d)",
                                                      px,py,batch);
                    }
                    return UAMMD_GET_3D_ROW_MAJOR(xy_height,Nbatches,Npx,Npy,batch,px,py);
                }

        };
    }

    template<class AFMimageType>
    class AFMimageInteractor: public Interactor{

        private:

            std::shared_ptr<GlobalData> gd;

            const int THREADS_PER_X = 8;
            const int THREADS_PER_Y = 8;
            const int THREADS_PER_N = 8;

            int Nparticles; // The total number of particles in the simulation

            int Npx;
            int Npy;
            int Nbatches;

            real sigma;
            real gamma;
            real tipRadius;

            std::shared_ptr<AFMimageType> potential;

            // The following data structures are used to store the information about the
            // set of reference images provided. One image per batch.

            std::map<int,real2> batchOffset;

            std::map<int,real2> batchMin;
            std::map<int,real2> batchMax;

            std::map<int,real2> batchResolution;

            std::map<int,std::map<std::pair<int,int>,real3>> batchId2pixels;


            // GPU data
            int imageSize; // Npx*Npy
            int imageBufferSize; // Npx*Npy*Nbatches
            thrust::host_vector<real> referenceImages_h;    // It stores the reference images for each batch, it has size (Npx*Npy*Nbatches)
            thrust::device_vector<real> referenceImages_d;  // Fro the reference images we have two vectors, one for the host and one for the device
                                                           // The host vector is used to store the reference images in CPU memory, it is useful
                                                           // since the reference images are constant and are known before the simulation starts (from input data)
            thrust::device_vector<real> generatedImages_d; // It stores the generated images for each batch, it has size (Npx*Npy*Nbatches)

            thrust::host_vector<real2> pixelPos_h; // It stores the position of each pixel, it has size (Npx*Npy*Nbatches)
            thrust::device_vector<real2> pixelPos_d; // It stores the position of each pixel, it has size (Npx*Npy*Nbatches)

            // All the previus buffers are plain arrays, we use UAMMD_SET_3D_ROW_MAJOR to set the values
            // and UAMMD_GET_3D_ROW_MAJOR to get the values.


            int bufferSize; // Npx*Npy*N
            thrust::device_vector<real4> buffer_d;
            thrust::device_vector<real4> bufferAux_d; // This buffer is used to store some values during the updateBuffer kernel
                                                      // this values do not fit in the buffer_d (it is just because we need more space)
            // First part: buffer stores the values of the phi function and its derivative
            // Second part: buffer stores the values of GEN, REF_GEN, REF_GENderivative, GEN_GENderivative

            // The size of the buffer is (N,Npx,Npy) and it stores information about individual particles
            // Note that here it uses the id of the particle, meaning that the buffer ignores the batch information itself

            //// Reduction indexing

            //                              Npx
            //                      ._____________________    _
            //                     /|                    /|    \
            //                    / |                   / |    |
            //                   /  |                  /  |    |
            //          Npy     /   |                 /   |    |
            //                 /    |                /    |    |
            //                /     |               /     |    |
            //               /_     |              /      |    |
            //              /_/___________________/       |    |
            //              | |     |             |       |    |
            //              | |     |             |       |    |
            //              | |     |             |       |    |
            //              | |     |_____________|_______|    |
            //              | |    /|             |      /|    |
            //  Nbatch_1    | |   / |             |     / |    |
            //              | |  /  |             |    /  |    |
            //              | | /   .             |   /   .    |
            //              | |/    .             |  /    .    |
            //              | |     .             | /     .    |
            //              |/|                   |/           |
            //              |_/___________________/            |
            //              |                     |            |
            //              |                     |            | N particles
            //              |                     |            |
            //              .                     .            |
            //              .                     .            |
            //              .       .             .       .    |
            //                      .                     .    |
            //                      .                     .    |
            //  Nbatch_n            |                     |    |
            //                      |                     |    |
            //                      |_____________________|    |
            //                     /|                    /|    |
            //              .     / |             .     / |    |
            //              .    /  |             .    /  |    |
            //              .   /   |             .   /   |    |
            //              |  /    |             |  /    |    |
            //              | /     |_____________|_/_____|    |
            //              |/     /|_____________|/____//|    |
            //              |_____________________/    // |    |
            //              |    // |             |   //  |    |
            //              |   //  |             |  //   |    |
            //              |  //   |             | //    |    |
            //              | //    |_____________|//_____|   _/
            //              |/____________________//     /
            //              |/____________________/     /
            //              |    /                |    /
            //              |   /                 |   /
            //              |  /                  |  /
            //              | /                   | /
            //              |/                    |/
            //              |_____________________/
            //

            // The following arrays are used to store the indices of the buffer are used in the reductions
            // There are two types of reductions:
            // 1. Reduction of the columns of the each batch of the buffer (see first cube in the figure)
            //    The result of this reduction is a 3D array of size (Nbatches,Npx,Npy)
            //    We call this reduction PixelReduction
            // 2. Reduction of the planes of the buffer (see final cube in the figure)
            //    The result of this reduction is a 1D array of size (Nparticles)
            //    We call this reduction ParticleImageReduction

            // The following arrays store the indices of the buffer that correspond to each pixel of the PixelReduction
            // and the starting index of each segment (each pixel) in the reduction.
            thrust::device_vector<int> pixelReductionIndices_d;
            thrust::device_vector<int> pixelReductionIndicesSegmentStart_d;

            // The following arrays store the indices of the buffer that correspond to each particle image of the ParticleImageReduction
            // and the starting index of each segment (each particle image) in the reduction.
            thrust::device_vector<int> particleImageReductionIndices_d;
            thrust::device_vector<int> particleImageReductionIndicesSegmentStart_d;

            // Thrust iterators for the reduction indices
            using ElementIterator = thrust::device_vector<real4>::iterator;
            using IndexIterator   = thrust::device_vector<int>::iterator;

            thrust::permutation_iterator<ElementIterator,IndexIterator> pixelReductionIterator;
            thrust::permutation_iterator<ElementIterator,IndexIterator> particleImageReductionIterator;
            thrust::permutation_iterator<ElementIterator,IndexIterator> particleImageReductionIteratorAux;

            void initializeReductionIterators(){
                pixelReductionIterator         = thrust::make_permutation_iterator<ElementIterator,IndexIterator>
                                                 (buffer_d.begin(),pixelReductionIndices_d.begin());

                particleImageReductionIterator = thrust::make_permutation_iterator<ElementIterator,IndexIterator>
                                                 (buffer_d.begin(),particleImageReductionIndices_d.begin());

                particleImageReductionIteratorAux = thrust::make_permutation_iterator<ElementIterator,IndexIterator>
                                                    (bufferAux_d.begin(),particleImageReductionIndices_d.begin());
            }

            // Cub temporary storage
            size_t cub_temp_storage_bytes;
            thrust::device_vector<std::uint8_t> cub_temp_storage_d; //uint8_t is an unsigned char

            ///// Reduction results

            // The following arrays store the results of the PixelReduction and ParticleImageReduction
            thrust::device_vector<real4> pixelReductionResult_d;
            thrust::device_vector<real4> particleImageReductionResult_d;
            thrust::device_vector<real4> particleImageReductionResultAux_d;

            // Auxiliary function to check the consistency of the reduction indices

            void CheckReductionIndices(const thrust::host_vector<int>& indices,
                                       const thrust::host_vector<int>& segmentStart,
                                       int Nsegments, // Number of segments
                                       int N, // Size of the buffer
                                       std::string name, std::string type){

                // Check if all elements of the indices have been filled
                // and that the last element is equal to the size of the buffer
                for(int i=0;i<indices.size();i++){
                    if(indices[i] == -1){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Error in the reduction indices (%s)",
                                                      name.c_str(),type.c_str());
                    }
                }

                for(int i=0;i<segmentStart.size();i++){
                    if(segmentStart[i] == -1){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Error in the reduction indices segment start,"
                                                      " some elements are not filled (%s), index=%d",
                                                      name.c_str(),type.c_str(),i);
                    }
                }

                if(segmentStart[Nsegments] != N){
                    System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Error in the reduction indices segment start,"
                                                  " the last element is not equal to the size of the buffer (%s)",
                                                  name.c_str(),type.c_str());
                }

                for(int i=0;i<segmentStart.size()-1;i++){
                    if(segmentStart[i] >= segmentStart[i+1]){
                        System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Error in the reduction indices segment start,"
                                                      " the elements are not in increasing order (%s)",
                                                      name.c_str(),type.c_str());
                    }
                }

            }

            // Auxiliary function to compute the size of the temporary storage needed by cub
            template<typename InputIteratorT, typename OutputIteratorT, typename OffsetIteratorT>
            size_t computeCubTempStorage(InputIteratorT d_input,
                                         OutputIteratorT d_output,
                                         int nSegments,
                                         OffsetIteratorT d_offsets){

                AFMimageInteractor_ns::Real4Sum reduction_op;
                real4 zero = make_real4(0.0);

                void* temp_storage_ptr = nullptr;
                size_t temp_storage_bytes = 0;

                cub::DeviceSegmentedReduce::Reduce(temp_storage_ptr,temp_storage_bytes,
                                                   d_input,d_output,nSegments,
                                                   d_offsets,d_offsets+1,
                                                   reduction_op,zero,0);
                cudaDeviceSynchronize();

                return temp_storage_bytes;
            }

            // Auxiliary function to perform the reduction
            template<typename InputIteratorT, typename OutputIteratorT, typename OffsetIteratorT>
            void performCubReduction(InputIteratorT d_input,
                                     OutputIteratorT d_output,
                                     int nSegments,
                                     OffsetIteratorT d_offsets,
                                     cudaStream_t st){

                AFMimageInteractor_ns::Real4Sum reduction_op;
                real4 zero = make_real4(0.0);

                void* temp_storage_ptr = thrust::raw_pointer_cast(cub_temp_storage_d.data());

                cub::DeviceSegmentedReduce::Reduce(temp_storage_ptr,cub_temp_storage_bytes,
                                                   d_input,d_output,nSegments,
                                                   d_offsets,d_offsets+1,
                                                   reduction_op,zero,st);
                cudaDeviceSynchronize();
            }


            void initializeCubTempStorage(){

                size_t cub_temp_storage_bytes_pixelReduction =
                computeCubTempStorage(pixelReductionIterator,
                                      thrust::raw_pointer_cast(pixelReductionResult_d.data()),
                                      imageBufferSize,
                                      thrust::raw_pointer_cast(pixelReductionIndicesSegmentStart_d.data()));

                size_t cub_temp_storage_bytes_particleImageReduction =
                computeCubTempStorage(particleImageReductionIterator,
                                      thrust::raw_pointer_cast(particleImageReductionResult_d.data()),
                                      Nparticles,
                                      thrust::raw_pointer_cast(particleImageReductionIndicesSegmentStart_d.data()));

                cub_temp_storage_bytes = std::max(cub_temp_storage_bytes_pixelReduction,
                                                  cub_temp_storage_bytes_particleImageReduction);

                cub_temp_storage_d.resize(cub_temp_storage_bytes);

                System::log<System::MESSAGE>("[AFMimageInteractor] (%s) Temporary storage for cub has been initialized (size=%d bytes)",
                                             name.c_str(),cub_temp_storage_bytes);
            }

            void performPixelReduction(cudaStream_t st){
                performCubReduction(pixelReductionIterator,
                                    thrust::raw_pointer_cast(pixelReductionResult_d.data()),
                                    imageBufferSize,
                                    thrust::raw_pointer_cast(pixelReductionIndicesSegmentStart_d.data()),
                                    st);
            }

            void performParticleImageReduction(cudaStream_t st){
                performCubReduction(particleImageReductionIterator,
                                    thrust::raw_pointer_cast(particleImageReductionResult_d.data()),
                                    Nparticles,
                                    thrust::raw_pointer_cast(particleImageReductionIndicesSegmentStart_d.data()),
                                    st);

                // We also need to perform the reduction for the bufferAux_d
                performCubReduction(particleImageReductionIteratorAux,
                                    thrust::raw_pointer_cast(particleImageReductionResultAux_d.data()),
                                    Nparticles,
                                    thrust::raw_pointer_cast(particleImageReductionIndicesSegmentStart_d.data()),
                                    st);
            }

            // Buffer related functions

            void populateBuffer(cudaStream_t st){

                // Populate the buffer
                dim3 threadsPerBlock(THREADS_PER_X,THREADS_PER_Y,THREADS_PER_N);
                dim3 numBlocks;

                numBlocks.x = Npx/threadsPerBlock.x + ((Npx%threadsPerBlock.x)?1:0);
                numBlocks.y = Npy/threadsPerBlock.y + ((Npy%threadsPerBlock.y)?1:0);
                numBlocks.z = Nparticles/threadsPerBlock.z + ((Nparticles%threadsPerBlock.z)?1:0);

                int* ids = this->pd->getId(access::location::gpu, access::mode::read).raw();

                real4* pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite).raw();
                int* batchId = this->pd->getBatchId(access::location::gpu, access::mode::readwrite).raw();

                real4* buffer = thrust::raw_pointer_cast(buffer_d.data());

                real2* pixelPos   = thrust::raw_pointer_cast(pixelPos_d.data());

                AFMimageInteractor_ns::populateBuffer<<<numBlocks,threadsPerBlock,0,st>>>(buffer,
                                                                                          ids,
                                                                                          pos,batchId,
                                                                                          sigma,gamma,tipRadius,
                                                                                          pixelPos,
                                                                                          Npx,Npy,Nparticles,Nbatches,
                                                                                          gd->getEnsemble()->getBox());
            }

            void updateBuffer(cudaStream_t st){
                dim3 threadsPerBlock(THREADS_PER_X,THREADS_PER_Y,THREADS_PER_N);
                dim3 numBlocks;

                numBlocks.x = Npx/threadsPerBlock.x + ((Npx%threadsPerBlock.x)?1:0);
                numBlocks.y = Npy/threadsPerBlock.y + ((Npy%threadsPerBlock.y)?1:0);
                numBlocks.z = Nparticles/threadsPerBlock.z + ((Nparticles%threadsPerBlock.z)?1:0);

                int* ids = this->pd->getId(access::location::gpu, access::mode::read).raw();

                real4* buffer = thrust::raw_pointer_cast(buffer_d.data());
                real4* bufferAux = thrust::raw_pointer_cast(bufferAux_d.data());
                real4* pixelReductionResult = thrust::raw_pointer_cast(pixelReductionResult_d.data());

                real*  referenceImages = thrust::raw_pointer_cast(referenceImages_d.data());

                int* batchId = this->pd->getBatchId(access::location::gpu, access::mode::readwrite).raw();

                //We call updateBuffer using current AFMimageType as template parameter
                AFMimageInteractor_ns::updateBuffer<AFMimageType><<<numBlocks,threadsPerBlock,0,st>>>(buffer,
                                                                                                      bufferAux,
                                                                                                      ids,
                                                                                                      pixelReductionResult,
                                                                                                      referenceImages,
                                                                                                      batchId,
                                                                                                      gamma,
                                                                                                      Npx,Npy,Nparticles,Nbatches);

            }

            void computeGeneratedImagesFromPixelReductionResults(cudaStream_t st){
                dim3 threadsPerBlock(THREADS_PER_X,THREADS_PER_Y,THREADS_PER_N);
                dim3 numBlocks;

                numBlocks.x = Npx/threadsPerBlock.x + ((Npx%threadsPerBlock.x)?1:0);
                numBlocks.y = Npy/threadsPerBlock.y + ((Npy%threadsPerBlock.y)?1:0);
                numBlocks.z = Nbatches/threadsPerBlock.z + ((Nbatches%threadsPerBlock.z)?1:0);

                real4* pixelReductionResult_ptr = thrust::raw_pointer_cast(pixelReductionResult_d.data());
                real*  generatedImages_ptr      = thrust::raw_pointer_cast(generatedImages_d.data());

                AFMimageInteractor_ns::computeGeneratedImagesFromPixelReductionResults<<<numBlocks,threadsPerBlock,0,st>>>(pixelReductionResult_ptr,
                                                                                                                           generatedImages_ptr,
                                                                                                                           gamma,
                                                                                                                           Npx,Npy,Nbatches);
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

                Nparticles = pg->getParticleData()->getNumParticles();
                if(Nparticles != pg->getNumberParticles()){
                    System::log<System::CRITICAL>("[AFMimageInteractor] (%s) The current implementation of AFMimageInteractor only supports"
                                                  " ParticleGroup that contains all the particles in the simulation",
                                                  name.c_str());
                }

                potential = std::make_shared<AFMimageType>(gd,pg,data);

                ////////////////// START OF INPUT DATA PREPARATION //////////////////

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

                //Set up all sizes
                imageSize = Npx*Npy;
                imageBufferSize = Npx*Npy*Nbatches;
                bufferSize = Npx*Npy*Nparticles;

                // Since we are assuming that each batch has the same number of pixels, Npx*Npy,
                // we can check that the number of pixels is consistent with the number of batches

                if(px.size() != imageBufferSize || py.size() != imageBufferSize || height.size() != imageBufferSize){
                    System::log<System::CRITICAL>("[AFMimageInteractor] (%s) Number of pixels is inconsistent with the given image size (Npx=%d, Npy=%d, total=%d) "
                                                  "and the number of batches (Nbatches=%d). The total given data should be Npx*Npy*Nbatches=%d",
                                                  name.c_str(),Npx,Npy,imageSize,Nbatches,imageBufferSize);
                }

                // At this point batchIds is consistent with the batch information in the ParticleData
                // and it is loaded in the vector batchIds


                // Check pixels
                // We create a data structure that associates each batch with another map that associates
                // each pixel (using int coordinates) with the pixel value (real3, x,y,height)
                // std::map<int,std::map<std::pair<int,int>,real3>> batchId2pixels; (Defined above)

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

                // At this point we have set:
                // batchOffset: where the image starts. But the position of the pixel is the center of the pixel
                // batchMin: the minimum position of the image. The true minimum position, batchOffset - resolution/2
                // batchMax: the maximum position of the image. The true maximum position, batchOffset + resolution/2 + resolution*(Npx-1,Npy-1)
                // batchResolution: the resolution of the image
                // batchId2pixels: a data structure that associates each batch with another map that associates
                //                each pixel (using int coordinates) with the pixel value (real3, x,y,height)

                // Now we can fill the referenceImages_h vector and the pixelPos_h vector
                referenceImages_h.resize(imageBufferSize);
                pixelPos_h.resize(imageBufferSize);
                for(int currentBatch=0;currentBatch<Nbatches;currentBatch++){
                    for(int i=0;i<Npx;i++){
                    for(int j=0;j<Npy;j++){
                        real3 pixel = batchId2pixels.at(currentBatch).at(std::make_pair(i,j)); // We use .at() to ensure that the key exists

                        UAMMD_SET_3D_ROW_MAJOR(referenceImages_h,Nbatches,Npx,Npy,currentBatch,i,j,pixel.z);
                        UAMMD_SET_3D_ROW_MAJOR(pixelPos_h,Nbatches,Npx,Npy,currentBatch,i,j,make_real2(pixel.x,pixel.y));
                    }}
                }

                ////////////////// END OF INPUT DATA PREPARATION ////////////////// (CPU)


                ////////////////// START OF GPU DATA PREPARATION //////////////////

                // Resize the device vectors
                referenceImages_d = referenceImages_h; // This vector is constant, we can copy it here
                generatedImages_d.resize(imageBufferSize);

                // Copy the pixelPos_h to the device
                pixelPos_d = pixelPos_h;

                // Resize the buffer_d
                buffer_d.resize(bufferSize);
                bufferAux_d.resize(bufferSize);

                // Resize the pixelReductionResult_d and particleImageReductionResult_d
                pixelReductionResult_d.resize(imageBufferSize);
                particleImageReductionResult_d.resize(Nparticles);
                particleImageReductionResultAux_d.resize(Nparticles);

                /////// Reduction indexing
                //// Prepare the pixelReductionIndices_d and pixelReductionIndicesSegmentStart_d

                thrust::host_vector<int> pixelReductionIndices_h(bufferSize,-1);
                thrust::host_vector<int> pixelReductionIndicesSegmentStart_h(imageBufferSize+1,-1);

                // Now we are going to create an auxiliary array that associates each batch with the
                // ids of the particles that belong to that batch.

                std::map<int,std::vector<int>> batchId2ids;
                for(int i=0;i<Nbatches;i++){
                    batchId2ids[i] = std::vector<int>();
                }

                auto ids     = this->pd->getId(access::location::cpu, access::mode::read);
                auto batchId = this->pd->getBatchId(access::location::cpu, access::mode::read);

                for(int i=0;i<Nparticles;i++){
                    batchId2ids[batchId[i]].push_back(ids[i]);
                }

                int pixelReductionIndex = 0;
                int currentSegmentStart = 0;

                // First we iterate (IN ORDER) over the elements of generatedImages_d
                for(int generatedImage_index=0;
                        generatedImage_index<imageBufferSize;
                        generatedImage_index++){

                    int3 coordinates = MatrixIndexing::rowMajor3Dcoordinates(Nbatches,Npx,Npy,generatedImage_index);

                    int currentBatch = coordinates.x;
                    int px = coordinates.y;
                    int py = coordinates.z;

                    // Now we have to determine the indices in the buffer that correspond to the pixel (currentBatch,px,py)
                    for(int i=0;i<batchId2ids[currentBatch].size();i++){
                        int id = batchId2ids[currentBatch][i];
                        int index = MatrixIndexing::rowMajor3Dindex(Nparticles,Npx,Npy,id,px,py);
                        pixelReductionIndices_h[pixelReductionIndex] = index;
                        pixelReductionIndex++;
                    }
                    pixelReductionIndicesSegmentStart_h[generatedImage_index] = currentSegmentStart;
                    currentSegmentStart += batchId2ids[currentBatch].size();
                }
                pixelReductionIndicesSegmentStart_h[imageBufferSize] = currentSegmentStart;

                CheckReductionIndices(pixelReductionIndices_h,
                                      pixelReductionIndicesSegmentStart_h,
                                      imageBufferSize,
                                      bufferSize,
                                      name,"PixelReduction");

                // Copy the data to the device
                pixelReductionIndices_d = pixelReductionIndices_h;
                pixelReductionIndicesSegmentStart_d = pixelReductionIndicesSegmentStart_h;

                //// Prepare the particleImageReductionIndices_d and particleImageReductionIndicesSegmentStart_d

                thrust::host_vector<int> particleImageReductionIndices_h(bufferSize,-1);
                thrust::host_vector<int> particleImageReductionIndicesSegmentStart_h(Nparticles+1,-1);

                int particleImageReductionIndex = 0;
                currentSegmentStart = 0;

                // First we iterate (IN ORDER) over the elements the ids of the particles
                for(int i=0;i<Nparticles;i++){ // We assume that i is the id of the particle

                    // Now we have to determine the indices in the buffer that correspond to the particle id
                    int currentBatch = batchId[i];
                    for(int px=0;px<Npx;px++){
                    for(int py=0;py<Npy;py++){
                        int index = MatrixIndexing::rowMajor3Dindex(Nparticles,Npx,Npy,i,px,py);
                        particleImageReductionIndices_h[particleImageReductionIndex] = index;
                        particleImageReductionIndex++;
                    }}
                    particleImageReductionIndicesSegmentStart_h[i] = currentSegmentStart;
                    currentSegmentStart += imageSize;
                }
                particleImageReductionIndicesSegmentStart_h[Nparticles] = currentSegmentStart;

                CheckReductionIndices(particleImageReductionIndices_h,
                                      particleImageReductionIndicesSegmentStart_h,
                                      Nparticles,
                                      bufferSize,
                                      name,"ParticleImageReduction");

                // Copy the data to the device
                particleImageReductionIndices_d = particleImageReductionIndices_h;
                particleImageReductionIndicesSegmentStart_d = particleImageReductionIndicesSegmentStart_h;

                ////
                // Initialize thrust iterators
                initializeReductionIterators();

                //// Initialize the cub temporary storage
                initializeCubTempStorage();

                ////////////////// END OF GPU DATA PREPARATION //////////////////
            }

            std::shared_ptr<AFMimageInteractor_ns::AFMimage> takeImage(cudaStream_t st){
                System::log<System::DEBUG>("[AFMimageInteractor] (%s) Taking image",
                                           name.c_str());
                {
                    // Populate the buffer
                    populateBuffer(st);
                    // Perform the pixel reduction
                    performPixelReduction(st);
                    // Compute the generated images from the pixel reduction results
                    computeGeneratedImagesFromPixelReductionResults(st);

                    cudaDeviceSynchronize();
                }

                thrust::host_vector<real> generatedImages_h = generatedImages_d;

                return std::make_shared<AFMimageInteractor_ns::AFMimage>(pixelPos_h,generatedImages_h,
                                                                         Npx,Npy,Nbatches);
            }

            void sum(Computables comp,cudaStream_t st) override {

                if(comp.energy == true || comp.force == true){
                    // Common for both energy and force
                    System::log<System::DEBUG>("[AFMimageInteractor] (%s) Computing common data",
                                               name.c_str());

                    // Populate the buffer
                    populateBuffer(st);

                    // Perform the pixel reduction
                    performPixelReduction(st);

                    // At this point we have the buffer with the values of the phi function and its derivative
                    // and we have the quantity sum_(i=1)^N_i phi_i for each pixel stored in the component w of pixelReductionResult_d
                    // (also the sum of the derivatives of the phi function, which is not used)

                    // FOR DEBUGGING:
                    //// Compute the generated images from the pixel reduction results
                    //computeGeneratedImagesFromPixelReductionResults(st);

                    //// At this point we have the generated images for each batch stored in the generatedImages_d vector
                    //// Note this is a result itself, but it is used to compute the force and energy
                    //// It can be checked using generatedImages_h = generatedImages_d; (for debugging)
                    // generatedImages_h = generatedImages_d;
                    // For pixel i,j of batch n we have to use UAMMD_GET_3D_ROW_MAJOR(generatedImages_h,Nbatches,Npx,Npy,n,i,j)
                    // END FOR DEBUGGING

                    // Update the buffer
                    updateBuffer(st);

                    // Perform the particle image reduction
                    performParticleImageReduction(st);
                    // The results of the reduction are stored in particleImageReductionResult_d and particleImageReductionResultAux_d
                    // That have size Nparticles
                }

                // We use one thread per particle
                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<AFMimageType>::value){

                        int Nthreads = 256;
                        int Nblocks  = Nparticles/Nthreads + ((Nparticles%Nthreads)?1:0);

                        int* ids = this->pd->getId(access::location::gpu, access::mode::readwrite).raw();

                        real4* particleImageReductionResult_ptr = thrust::raw_pointer_cast(particleImageReductionResult_d.data());
                        real4* particleImageReductionResultAux_ptr = thrust::raw_pointer_cast(particleImageReductionResultAux_d.data());

                        AFMimageInteractor_ns::transverseAFMimage
                        <typename AFMimageType::ComputationalData,
                         typename AFMimageType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(potential->getComputationalData(comp,st),
                                                     potential->getEnergyTransverser(),
                                                     ids,
                                                     particleImageReductionResult_ptr,
                                                     particleImageReductionResultAux_ptr,
                                                     Nparticles);

                        CudaCheckError();


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

                        int Nthreads = 256;
                        int Nblocks  = Nparticles/Nthreads + ((Nparticles%Nthreads)?1:0);

                        int* ids = this->pd->getId(access::location::gpu, access::mode::readwrite).raw();

                        real4* particleImageReductionResult_ptr = thrust::raw_pointer_cast(particleImageReductionResult_d.data());
                        real4* particleImageReductionResultAux_ptr = thrust::raw_pointer_cast(particleImageReductionResultAux_d.data());

                        AFMimageInteractor_ns::transverseAFMimage
                        <typename AFMimageType::ComputationalData,
                         typename AFMimageType::ForceTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(potential->getComputationalData(comp,st),
                                                     potential->getForceTransverser(),
                                                     ids,
                                                     particleImageReductionResult_ptr,
                                                     particleImageReductionResultAux_ptr,
                                                     Nparticles);

                        CudaCheckError();

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
