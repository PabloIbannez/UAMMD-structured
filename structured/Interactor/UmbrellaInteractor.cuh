#ifndef __UMBRELLA_INTERACTOR__
#define __UMBRELLA_INTERACTOR__

namespace uammd{
namespace structured{
namespace Interactor{

    namespace UmbrellaSetInteractor_ns{

        template<class umbInfoType,int BLOCK_DIM>
        void __global__ umbrellaSetInteractorForceKernel(real*  mass,
                                                         real4* pos,
                                                         real4* force,
                                                         umbInfoType* ui,
                                                         int Nsets,
                                                         Box box,
                                                         real K){
            
            typedef cub::BlockReduce<real3,BLOCK_DIM> BlockReduce;
            
            int set = blockIdx.x;
            if(set >= Nsets){return;}
            
            __shared__ int N1;
            __shared__ int N2;
            __shared__ const int* g1;
            __shared__ const int* g2;
            __shared__ real M1;
            __shared__ real M2;
            __shared__ real r0;
            
            __shared__ real3 F12;

            if(threadIdx.x == 0){
                N1 = ui[set].N1;
                N2 = ui[set].N2;
                g1 = ui[set].g1;
                g2 = ui[set].g2;
                M1 = ui[set].M1;
                M2 = ui[set].M2;
                r0 = ui[set].r0;
            }
            __syncthreads();

            __shared__ typename BlockReduce::TempStorage temp_storage;

            if(threadIdx.x < N1 or threadIdx.x < N2){
                
                real3 m1p1=make_real3(0.0);
                for(int i=threadIdx.x;i<N1;i+=blockDim.x){
                    m1p1+=make_real3(pos[g1[i]])*mass[g1[i]];
                }
                
                real3 m2p2=make_real3(0.0);
                for(int i=threadIdx.x;i<N2;i+=blockDim.x){
                    m2p2+=make_real3(pos[g2[i]])*mass[g2[i]];
                }
            
                real3 COM1 = BlockReduce(temp_storage).Sum(m1p1,N1)/M1;
                __syncthreads();
                
                real3 COM2 = BlockReduce(temp_storage).Sum(m2p2,N2)/M2;
            
                if(threadIdx.x == 0){
                    const real3 r12 = box.apply_pbc(COM2-COM1);

                    const real  r2  = dot(r12,r12);
                    const real invr = rsqrt(r2);
                
                    //Force over 1
                    F12 =  K*(real(1.0)-r0*invr)*r12;
                }
                __syncthreads();
                
                //Distr forces
                for(int i=threadIdx.x;i<N1;i+=blockDim.x){
                    force[g1[i]] += make_real4(F12*mass[g1[i]]/M1,real(0.0));
                }
                for(int i=threadIdx.x;i<N2;i+=blockDim.x){
                    force[g2[i]] += make_real4(-F12*mass[g2[i]]/M2,real(0.0));
                }
            }
        }
        
        template<class umbInfoType>
        void __global__ umbrellaSetInteractorEnergyKernel(real*  mass,
                                                          real4* pos,
                                                          real*  energy,
                                                          umbInfoType* ui,
                                                          int Nsets,
                                                          Box box,
                                                          real K){}
    }

    
    template<int BLOCK_DIM=128>
    class UmbrellaSetInteractor: public Interactor{

        public:
            
            struct umbrellaInfo{
                int N1;
                int N2;
                const int* g1;
                const int* g2;
                real M1;
                real M2;
                real r0;
            };

        protected:
                
            int Nsets;

            real K;

            Box box;
                                  
            std::vector<std::shared_ptr<ParticleGroup>> pg1set;
            std::vector<std::shared_ptr<ParticleGroup>> pg2set;

            thrust::device_vector<umbrellaInfo> umbInfo;
            thrust::host_vector<umbrellaInfo> umbInfo_h;

            connection reorderConnection;

            void updateUmbrellaInfo(){
                fori(0,Nsets){
                    umbInfo_h[i].g1 = pg1set[i]->getIndicesRawPtr(access::location::gpu);
                    umbInfo_h[i].g2 = pg2set[i]->getIndicesRawPtr(access::location::gpu);
                }
                umbInfo = umbInfo_h;
            }

        public:

            struct Parameters{
                real K;
            };
            
            UmbrellaSetInteractor(std::shared_ptr<ParticleGroup> pg,
                                  std::vector<std::shared_ptr<ParticleGroup>> pg1set,
                                  std::vector<std::shared_ptr<ParticleGroup>> pg2set,
                                  std::vector<real> r0,
                                  Parameters par):Interactor(pg,std::string("UmbrellaSetInteractor")),
                                                  K(par.K),
                                                  pg1set(pg1set),
                                                  pg2set(pg2set){
              
                reorderConnection  = pd->getReorderSignal()->connect([this](){this->updateUmbrellaInfo();});

                if(!(pg1set.size() == pg2set.size()) or !(pg1set.size() == r0.size())){
                    this->sys->template log<uammd::System::CRITICAL>("[UmbrellaSetInteractor] The size of umbrella sets,and the r0 vector "
                                                                      "have to be equal size, but sizes are %i %i (sets) %i (r0).",
                                                                      pg1set.size(),pg2set.size(),r0.size());
                } else {
                    Nsets = pg1set.size();
                }

                umbInfo.resize(Nsets);
                umbInfo_h.resize(Nsets);

                fori(0,Nsets){
                    umbInfo_h[i].N1 = pg1set[i]->getNumberParticles();
                    umbInfo_h[i].N2 = pg2set[i]->getNumberParticles();
                    
                    umbInfo_h[i].M1 = Measures::totalMass(sys,pd,pg1set[i],0);
                    umbInfo_h[i].M2 = Measures::totalMass(sys,pd,pg2set[i],0);

                    umbInfo_h[i].r0 = r0[i];
                }

                this->updateUmbrellaInfo();
            }

            ~UmbrellaSetInteractor(){
                  reorderConnection.disconnect();
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                
                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = BLOCK_DIM;
                int Nblocks  = Nsets;

                if(comp.force == true){
                    
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
                    
                    UmbrellaSetInteractor_ns::umbrellaSetInteractorForceKernel<umbrellaInfo,BLOCK_DIM><<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                                                     pos,
                                                                                                                                     force,
                                                                                                                                     thrust::raw_pointer_cast(umbInfo.data()),
                                                                                                                                     Nsets,
                                                                                                                                     box,
                                                                                                                                     K);  

                }
                
                if(comp.energy == true){
                    
                    auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
                    
                    UmbrellaSetInteractor_ns::umbrellaSetInteractorEnergyKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                              pos,
                                                                                                              energy,
                                                                                                              thrust::raw_pointer_cast(umbInfo.data()),
                                                                                                              Nsets,
                                                                                                              box,
                                                                                                              K);  

                }
                
                CudaCheckError();
            
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
    namespace UmbrellaAlongVectorSetInteractor_ns{

        template<class umbInfoType,int BLOCK_DIM>
        void __global__ umbrellaAlongVectorSetInteractorForceKernel(real*  mass,
                                                                    real4* pos,
                                                                    real4* force,
                                                                    umbInfoType* ui,
                                                                    int Nsets,
                                                                    Box box,
                                                                    real K){
            
            typedef cub::BlockReduce<real3,BLOCK_DIM> BlockReduce;
            
            int set = blockIdx.x;
            if(set >= Nsets){return;}
            
            __shared__ int N1;
            __shared__ const int* g1;
            __shared__ real M1;
            __shared__ real3 equiPos;
            
            __shared__ real3 F12;

            if(threadIdx.x == 0){
                N1 = ui[set].N1;
                g1 = ui[set].g1;
                M1 = ui[set].M1;
                equiPos = ui[set].equiPos;
            }
            __syncthreads();

            __shared__ typename BlockReduce::TempStorage temp_storage;

            if(threadIdx.x < N1){
                
                real3 m1p1=make_real3(0.0);
                for(int i=threadIdx.x;i<N1;i+=blockDim.x){
                    m1p1+=make_real3(pos[g1[i]])*mass[g1[i]];
                }
            
                real3 COM1 = BlockReduce(temp_storage).Sum(m1p1,N1)/M1;
            
                if(threadIdx.x == 0){
                    const real3 r12 = box.apply_pbc(equiPos-COM1);

                    const real  r2  = dot(r12,r12);
                    const real invr = rsqrt(r2);
                
                    //Force over 1
                    F12 =  K*r12;
                }
                __syncthreads();
                
                //Distr forces
                for(int i=threadIdx.x;i<N1;i+=blockDim.x){
                    force[g1[i]] += make_real4(F12*mass[g1[i]]/M1,real(0.0));
                }
            }
        }
        
        template<class umbInfoType>
        void __global__ umbrellaAlongVectorSetInteractorEnergyKernel(real*  mass,
                                                                     real4* pos,
                                                                     real*  energy,
                                                                     umbInfoType* ui,
                                                                     int Nsets,
                                                                     Box box,
                                                                     real K){}
    }

    
    template<int BLOCK_DIM=128>
    class UmbrellaAlongVectorSetInteractor: public Interactor{

        public:
            
            struct umbrellaInfo{
                int N1;
                const int* g1;
                real  M1;
                real3 equiPos;
            };

        protected:
                
            int Nsets;

            real K;

            Box box;
                                  
            std::vector<std::shared_ptr<ParticleGroup>> pg1set;

            thrust::device_vector<umbrellaInfo> umbInfo;
            thrust::host_vector<umbrellaInfo> umbInfo_h;

            connection reorderConnection;

            void updateUmbrellaInfo(){
                fori(0,Nsets){
                    umbInfo_h[i].g1 = pg1set[i]->getIndicesRawPtr(access::location::gpu);
                }
                umbInfo = umbInfo_h;
            }

        public:

            struct Parameters{
                real K;
            };
            
            UmbrellaAlongVectorSetInteractor(std::shared_ptr<ParticleGroup> pg,
                                             std::vector<std::shared_ptr<ParticleGroup>> pg1set,
                                             std::vector<real3> equiPos,
                                             Parameters par):Interactor(pg,std::string("UmbrellaSetInteractor")),
                                                             K(par.K),
                                                             pg1set(pg1set){
              
                reorderConnection  = pd->getReorderSignal()->connect([this](){this->updateUmbrellaInfo();});

                if(!(pg1set.size() == equiPos.size())){
                    this->sys->template log<uammd::System::CRITICAL>("[UmbrellaSetInteractor] The size of umbrella set,and the equiPos vector "
                                                                      "have to be equal size, but sizes are %i (set) %i (equiPos).",
                                                                      pg1set.size(),equiPos.size());
                } else {
                    Nsets = pg1set.size();
                }

                umbInfo.resize(Nsets);
                umbInfo_h.resize(Nsets);

                fori(0,Nsets){
                    umbInfo_h[i].N1 = pg1set[i]->getNumberParticles();
                    umbInfo_h[i].M1 = Measures::totalMass(sys,pd,pg1set[i],0);
                    umbInfo_h[i].equiPos = equiPos[i];
                }

                this->updateUmbrellaInfo();
            }

            ~UmbrellaAlongVectorSetInteractor(){
                  reorderConnection.disconnect();
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                
                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = BLOCK_DIM;
                int Nblocks  = Nsets;

                if(comp.force == true){
                    
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
                    
                    UmbrellaAlongVectorSetInteractor_ns::umbrellaAlongVectorSetInteractorForceKernel<umbrellaInfo,BLOCK_DIM><<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                                                                           pos,
                                                                                                                                                           force,
                                                                                                                                                           thrust::raw_pointer_cast(umbInfo.data()),
                                                                                                                                                           Nsets,
                                                                                                                                                           box,
                                                                                                                                                           K);  

                }
                
                if(comp.energy == true){
                    
                    auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
                    
                    UmbrellaAlongVectorSetInteractor_ns::umbrellaAlongVectorSetInteractorEnergyKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                                                    pos,
                                                                                                                                    energy,
                                                                                                                                    thrust::raw_pointer_cast(umbInfo.data()),
                                                                                                                                    Nsets,
                                                                                                                                    box,
                                                                                                                                    K);  

                }
                
                CudaCheckError();
            
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
}}}


#endif
