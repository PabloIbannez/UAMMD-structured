#ifndef __PLATES_INTERACTOR__
#define __PLATES_INTERACTOR__

namespace uammd{
namespace structured{

    namespace Interactor{
    namespace Plates_ns{

        void __global__ computeForceKernel(real4* pos,
                                           real4* force,
                                           real4* TopPlateForceTemp_ptr,
                                           real4* BottomPlateForceTemp_ptr,
                                           ParticleGroup::IndexIterator groupIter,
                                           real   TopPlatePos,
                                           real   BottomPlatePos,
                                           int numberParticles){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= numberParticles) return;

            const real3 p  = make_real3(pos[groupIter[i]]);

            const real dzTop    = TopPlatePos-p.z;
            
            const real invdzTop2  = real(1.0)/(dzTop*dzTop);
            const real invdzTop4  = invdzTop2*invdzTop2;
            const real invdzTop10 = invdzTop4*invdzTop4*invdzTop2;
                    
            const real3 fTop = {0,0,-invdzTop10*invdzTop2*dzTop};
            
            const real dzBottom    = BottomPlatePos-p.z;
            
            const real invdzBottom2  = real(1.0)/(dzBottom*dzBottom);
            const real invdzBottom4  = invdzBottom2*invdzBottom2;
            const real invdzBottom10 = invdzBottom4*invdzBottom4*invdzBottom2;
                    
            const real3 fBottom = {0,0,-invdzBottom10*invdzBottom2*dzBottom};

            const real4 fTotal = make_real4(fTop+fBottom,0);
            //Force on the particles due to the tip
            force[groupIter[i]] += fTotal;
            
            TopPlateForceTemp_ptr[i]    = make_real4(fTop,0);    //NOT +=
            BottomPlateForceTemp_ptr[i] = make_real4(fBottom,0); //NOT +=
        }
    }
    
    class Plates: public Interactor{

        private:

            int N;
                
            thrust::device_vector<real4> TopPlateForceTemp;
            thrust::device_vector<real4> BottomPlateForceTemp;
            
            real4* TopPlateTotalForce;
            real4* BottomPlateTotalForce;

                
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;
            
            real TopPlatePos;
            real BottomPlatePos;
            
        public:

            struct Parameters{};

            Parameters inputFileToParam(InputFile& in){

                Parameters param;

                return param;

            }

            Plates(std::shared_ptr<ParticleGroup> pg,
                   InputFile& in):Plates(pg,inputFileToParam(in)){}

            Plates(std::shared_ptr<ParticleGroup> pg,
                   Parameters par):Interactor(pg,std::string("Plates"))
            {
                
                N=pg->getNumberParticles();
                
                TopPlateForceTemp.resize(N);
                BottomPlateForceTemp.resize(N);

                cudaMallocManaged(&TopPlateTotalForce,    sizeof(real4));
                cudaMallocManaged(&BottomPlateTotalForce, sizeof(real4));

                cubReductionTempStorage     = NULL;
                cubReductionTempStorageSize = 0;
                
                cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                       cubReductionTempStorageSize,
                                       thrust::raw_pointer_cast(TopPlateForceTemp.data()), 
                                       TopPlateTotalForce, 
                                       N);

                cudaMalloc(&cubReductionTempStorage, cubReductionTempStorageSize);

                cudaDeviceSynchronize();
            }
            
            ~Plates(){
                cudaFree(cubReductionTempStorage);
                
                cudaFree(TopPlateTotalForce);
                cudaFree(BottomPlateTotalForce);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                    
                int blocksize = 256;
                int Nthreads = blocksize<N?blocksize:N;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                if(comp.force == true){
                
                    auto pos   = pd->getPos(access::location::gpu, access::mode::read).raw();
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    auto TopPlateForceTemp_ptr    = thrust::raw_pointer_cast(TopPlateForceTemp.data());
                    auto BottomPlateForceTemp_ptr = thrust::raw_pointer_cast(BottomPlateForceTemp.data());

                    auto groupIter = pg->getIndexIterator(access::location::gpu);
                                               
                    Plates_ns::computeForceKernel<<<Nblocks, Nthreads, 0, st>>>(pos,
                                                                                force,
                                                                                TopPlateForceTemp_ptr,
                                                                                BottomPlateForceTemp_ptr,
                                                                                groupIter,
                                                                                TopPlatePos,
                                                                                BottomPlatePos,
                                                                                N);  
                    CudaCheckError();

                }

                if(comp.energy == true){}
                if(comp.virial == true){}
                if(comp.stress == true){}

            }
            
            void    setTopPlatePosition(double newPos){TopPlatePos=newPos;}
            double  getTopPlatePosition(){return TopPlatePos;}
            
            void    setBottomPlatePosition(double newPos){BottomPlatePos=newPos;}
            double  getBottomPlatePosition(){return BottomPlatePos;}
            
            real2 getPlatesForce(cudaStream_t st){
                
                cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                       cubReductionTempStorageSize,
                                       thrust::raw_pointer_cast(TopPlateForceTemp.data()), 
                                       TopPlateTotalForce, 
                                       N,st);
                
                cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                       cubReductionTempStorageSize,
                                       thrust::raw_pointer_cast(BottomPlateForceTemp.data()), 
                                       BottomPlateTotalForce, 
                                       N,st);

                return {TopPlateTotalForce[0].z,BottomPlateTotalForce[0].z};
            }
    };

}}}

#endif
