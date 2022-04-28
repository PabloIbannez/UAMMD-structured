#ifndef __SETS_INTERACTOR__
#define __SETS_INTERACTOR__

namespace uammd{
namespace structured{
namespace Interactor{
namespace Sets{
    
    namespace ExternalForceOverCOM_ns{

        void __global__ computeExternalForceOverCOM_ForceKernel(const int * __restrict__ id2index,
                                                                real* mass,
                                                                real4* pos,
                                                                real4* force,
                                                                int  setSize,
                                                                int  nSets,
                                                                int* set2id,
                                                                real3* F){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            real totalMass = 0;
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];
                totalMass+=mass[index];
            }
 
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];

                const real3 f =  (mass[index]/totalMass)*F[set];
                force[index] += make_real4(f,0);
            }
        }
    }
    
    class ExternalForceOverCOM: public Interactor{

        protected:

            int nSets;
            int setSize;
            
            thrust::device_vector<int> set2id;
            
            thrust::device_vector<real3> externalForce;

        public:

            struct Parameters{};

        public:
            
            ExternalForceOverCOM(std::shared_ptr<System>       sys,
                                 std::shared_ptr<ParticleData>  pd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 int setSize,
                                 int nSets,
                                 thrust::host_vector<int>   set2id,
                                 thrust::host_vector<real3> externalForce,
                                 Parameters par):Interactor(pd,pg,sys,std::string("ExternalForceOverCOM")),
                                                 setSize(setSize),
                                                 nSets(nSets),
                                                 set2id(set2id),
                                                 externalForce(externalForce){}
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
    
                    ExternalForceOverCOM_ns::computeExternalForceOverCOM_ForceKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                                   mass,
                                                                                                                   pos,
                                                                                                                   force,
                                                                                                                   setSize,
                                                                                                                   nSets,
                                                                                                                   thrust::raw_pointer_cast(set2id.data()), 
                                                                                                                   thrust::raw_pointer_cast(externalForce.data()));  

                }
            }
    };
    
    namespace ExternalTorqueOverCOM_ns{

        void __global__ computeExternalTorqueOverCOM_ForceKernel(const int * __restrict__ id2index,
                                                                 real* mass,
                                                                 real4* pos,
                                                                 real4* force,
                                                                 int  setSize,
                                                                 int  nSets,
                                                                 int* set2id,
                                                                 real3* torque){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            //Compute center of mass
            real totalMass = real(0);
            real3 com = make_real3(0.0);
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];
                totalMass+=mass[index];
                com+=mass[index]*make_real3(pos[index]);
            }
            com=com/totalMass;
 
            //Axis
            real prefactor = real(0);
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];
                
                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com;
                
                real3 T  = torque[set];
                real3 ui = ri-T*dot(ri,T)/dot(T,T);
                prefactor+=mi*dot(ui,ui);
            }
            prefactor=real(1.0)/prefactor;
            
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];

                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com;
                
                real3 T  = torque[set];
                real3 ui = ri-T*dot(ri,T)/dot(T,T);
                
                const real3 f =  mi*prefactor*cross(T,ui);

                force[index] += make_real4(f,0);
            }
        }
    }
    
    class ExternalTorqueOverCOM: public Interactor{

        protected:

            int nSets;
            int setSize;
            
            thrust::device_vector<int> set2id;
            
            thrust::device_vector<real3> externalTorque;

        public:

            struct Parameters{};

        public:
            
            ExternalTorqueOverCOM(std::shared_ptr<System>       sys,
                                  std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<ParticleGroup> pg,
                                  int setSize,
                                  int nSets,
                                  thrust::host_vector<int>   set2id,
                                  thrust::host_vector<real3> externalTorque,
                                  Parameters par):Interactor(pd,pg,sys,std::string("ExternalTorqueOverCOM")),
                                                  setSize(setSize),
                                                  nSets(nSets),
                                                  set2id(set2id),
                                                  externalTorque(externalTorque){
                if(setSize == 1){
                    sys->log<uammd::System::CRITICAL>("[ExternalTorqueOverCOM] "
                                                      "The size of the group over the torque has to be applied is equal to one.");
                }
            }
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
    
                    ExternalTorqueOverCOM_ns::computeExternalTorqueOverCOM_ForceKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                                     mass,
                                                                                                                     pos,
                                                                                                                     force,
                                                                                                                     setSize,
                                                                                                                     nSets,
                                                                                                                     thrust::raw_pointer_cast(set2id.data()), 
                                                                                                                     thrust::raw_pointer_cast(externalTorque.data()));  
                }
            }
    };

    namespace HarmonicBondBtwCOM_ns{
        
        void __global__ computeHarmonicBondBtwCOM_ForceKernel(const int * __restrict__ id2index,
                                                               real* mass,
                                                               real4* pos,
                                                               real4* force,
                                                               int setSize1,
                                                               int setSize2,
                                                               int nSets,
                                                               int* set2id1,
                                                               int* set2id2,
                                                               real* r0,
                                                               real* K,
                                                               Box box){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            real totalMass1 = 0;
            real3 com1 = make_real3(0.0);

            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                totalMass1+=mass[index];
                com1+=make_real3(pos[index])*mass[index];
            }
            com1=com1/totalMass1;

            real totalMass2 = 0;
            real3 com2 = make_real3(0.0);

            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                totalMass2+=mass[index];
                com2+=make_real3(pos[index])*mass[index];
            }
            com2=com2/totalMass2;
 
            const real3 r12 = box.apply_pbc(com2-com1);
            
            const real  r2  = dot(r12,r12);
            const real invr = rsqrt(r2);
            
            //Force over 1
            real3 F12 =  Potentials::CommonPotentials::Harmonic::force(r12,r2,K[set],r0[set]);

            if(r0[set] == real(0.0) and r2 == real(0.0)){
                real3 F12 = make_real3(0.0);
            }
            
            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                const real3 f =  (mass[index]/totalMass1)*F12;
                force[index] += make_real4(f,0);
            }

            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                const real3 f =  -(mass[index]/totalMass2)*F12;
                force[index] += make_real4(f,0);
            }
        }
    }
    
    class HarmonicBondBtwCOM: public Interactor{

        protected:

            int nSets;
            
            int setSize1;
            int setSize2;
            
            thrust::device_vector<int> set2id1;
            thrust::device_vector<int> set2id2;
            
            thrust::device_vector<real> r0;
            thrust::device_vector<real> K;
            
            Box box;

        public:

            struct Parameters{};
            
        public:
            
            HarmonicBondBtwCOM(std::shared_ptr<System>       sys,
                                std::shared_ptr<ParticleData>  pd,
                                std::shared_ptr<ParticleGroup> pg,
                                int setSize1,
                                int setSize2,
                                int nSets,
                                thrust::host_vector<int> set2id1,
                                thrust::host_vector<int> set2id2,
                                thrust::host_vector<real> r0,
                                thrust::host_vector<real> K,
                                Parameters par):Interactor(pd,pg,sys,std::string("HarmonicBondBtwCOM")),
                                                setSize1(setSize1),
                                                setSize2(setSize2),
                                                nSets(nSets),
                                                set2id1(set2id1),
                                                set2id2(set2id2),
                                                r0(r0),
                                                K(K){}
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    HarmonicBondBtwCOM_ns::computeHarmonicBondBtwCOM_ForceKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                                 mass,
                                                                                                                 pos,
                                                                                                                 force,
                                                                                                                 setSize1, 
                                                                                                                 setSize2, 
                                                                                                                 nSets, 
                                                                                                                 thrust::raw_pointer_cast(set2id1.data()), 
                                                                                                                 thrust::raw_pointer_cast(set2id2.data()), 
                                                                                                                 thrust::raw_pointer_cast(r0.data()),  
                                                                                                                 thrust::raw_pointer_cast(K.data()),  
                                                                                                                 box);  

                }
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
    namespace HarmonicFixedCOM_ns{

        void __global__ computeHarmonicFixedCOM_ForceKernel(const int * __restrict__ id2index,
                                                            real* mass,
                                                            real4* pos,
                                                            real4* force,
                                                            int setSize,
                                                            int nSets,
                                                            int* set2id,
                                                            real3* fixedPoint,
                                                            real3* K,
                                                            Box box){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            real totalMass = 0;
            real3 com = make_real3(0.0);

            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];
                totalMass+=mass[index];
                com+=make_real3(pos[index])*mass[index];
            }
            com=com/totalMass;

            const real3 r1f = box.apply_pbc(fixedPoint[set]-com);

            //Force over 1
            real3 F1f =  Potentials::CommonPotentials::HarmonicAnisotropic::force(r1f,
                                                                                  K[set],
                                                                                  make_real3(0.0));
            for(int i=0;i<setSize;i++){
                int index = id2index[set2id[i+set*setSize]];
                const real3 f =  (mass[index]/totalMass)*F1f;
                force[index] += make_real4(f,0);
            }
        }
    }
    
    class HarmonicFixedCOM: public Interactor{

        protected:

            int nSets;
            int setSize;
            
            thrust::device_vector<int> set2id;
            
            thrust::device_vector<real3> fixedPoint;
            thrust::device_vector<real3> K;
                
            Box box;

        public:

            struct Parameters{};

        public:
            
            HarmonicFixedCOM(std::shared_ptr<System>       sys,
                             std::shared_ptr<ParticleData>  pd,
                             std::shared_ptr<ParticleGroup> pg,
                             int setSize,
                             int nSets,
                             thrust::host_vector<int> set2id,
                             thrust::host_vector<real3> fixedPoint,
                             thrust::host_vector<real3> K,
                             Parameters par):Interactor(pd,pg,sys,std::string("HarmonicFixedCOM")),
                                             setSize(setSize),
                                             nSets(nSets),
                                             set2id(set2id),
                                             fixedPoint(fixedPoint),
                                             K(K){}
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    HarmonicFixedCOM_ns::computeHarmonicFixedCOM_ForceKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                           mass,
                                                                                                           pos,
                                                                                                           force,
                                                                                                           setSize, 
                                                                                                           nSets, 
                                                                                                           thrust::raw_pointer_cast(set2id.data()), 
                                                                                                           thrust::raw_pointer_cast(fixedPoint.data()),  
                                                                                                           thrust::raw_pointer_cast(K.data()),  
                                                                                                           box);  
                }
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
    namespace ConstantForceBtwCOM_ns{
        
        void __global__ computeConstantForceBtwCOM_ForceKernel(const int * __restrict__ id2index,
                                                                real* mass,
                                                                real4* pos,
                                                                real4* force,
                                                                int setSize1,
                                                                int setSize2,
                                                                int nSets,
                                                                int* set2id1,
                                                                int* set2id2,
                                                                real* F,
                                                                Box box){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            real totalMass1 = 0;
            real3 com1 = make_real3(0.0);

            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                totalMass1+=mass[index];
                com1+=make_real3(pos[index])*mass[index];
            }
            com1=com1/totalMass1;

            real totalMass2 = 0;
            real3 com2 = make_real3(0.0);

            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                totalMass2+=mass[index];
                com2+=make_real3(pos[index])*mass[index];
            }
            com2=com2/totalMass2;
            
            const real3 r12 = box.apply_pbc(com2-com1);
            
            const real  r2  = dot(r12,r12);
            const real invr = rsqrt(r2);
            
            //Force over 1
            real3 F12 =  -F[set]*r12*invr;
            
            if(r2 == real(0.0)){
                real3 F12 = make_real3(0.0);
            }
            
            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                const real3 f =  (mass[index]/totalMass1)*F12;
                force[index] += make_real4(f,0);
            }

            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                const real3 f =  -(mass[index]/totalMass2)*F12;
                force[index] += make_real4(f,0);
            }
        }
    }
    
    class ConstantForceBtwCOM: public Interactor{

        protected:

            int nSets;
            
            int setSize1;
            int setSize2;
            
            thrust::device_vector<int> set2id1;
            thrust::device_vector<int> set2id2;
            
            thrust::device_vector<real> F;
            
            Box box;

        public:

            struct Parameters{};
            
            ConstantForceBtwCOM(std::shared_ptr<System>       sys,
                                 std::shared_ptr<ParticleData>  pd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 int setSize1,
                                 int setSize2,
                                 int nSets,
                                 thrust::host_vector<int> set2id1,
                                 thrust::host_vector<int> set2id2,
                                 thrust::host_vector<real> F,
                                 Parameters par):Interactor(pd,pg,sys,std::string("ConstantForceBtwCOM")),
                                                 setSize1(setSize1),
                                                 setSize2(setSize2),
                                                 nSets(nSets),
                                                 set2id1(set2id1),
                                                 set2id2(set2id2),
                                                 F(F){}
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    ConstantForceBtwCOM_ns::computeConstantForceBtwCOM_ForceKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                                   mass,
                                                                                                                   pos,
                                                                                                                   force,
                                                                                                                   setSize1, 
                                                                                                                   setSize2, 
                                                                                                                   nSets, 
                                                                                                                   thrust::raw_pointer_cast(set2id1.data()), 
                                                                                                                   thrust::raw_pointer_cast(set2id2.data()), 
                                                                                                                   thrust::raw_pointer_cast(F.data()),  
                                                                                                                   box);  

                }
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
    namespace ConstantTorqueBtwCOM_ns{
        
        void __global__ computeConstantTorqueBtwCOM_TorqueKernel(const int * __restrict__ id2index,
                                                                real* mass,
                                                                real4* pos,
                                                                real4* force,
                                                                int setSize1,
                                                                int setSize2,
                                                                int nSets,
                                                                int* set2id1,
                                                                int* set2id2,
                                                                real* torque,
                                                                Box box){

            int set = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(set >= nSets) return;

            real totalMass1 = 0;
            real3 com1 = make_real3(0.0);

            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                totalMass1+=mass[index];
                com1+=make_real3(pos[index])*mass[index];
            }
            com1=com1/totalMass1;

            real totalMass2 = 0;
            real3 com2 = make_real3(0.0);

            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                totalMass2+=mass[index];
                com2+=make_real3(pos[index])*mass[index];
            }
            com2=com2/totalMass2;
            
            real3 torqueVector = torque[set]*normalize(box.apply_pbc(com2-com1));
            
            ////////////////////////////

            //Axis
            real prefactor = real(0);
            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];
                
                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com1;
                
                real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);
                prefactor+=mi*dot(ui,ui);
            }
            prefactor=real(1.0)/prefactor;
            
            for(int i=0;i<setSize1;i++){
                int index = id2index[set2id1[i+set*setSize1]];

                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com1;
                
                real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);
                
                const real3 f =  mi*prefactor*cross(torqueVector,ui);

                force[index] += make_real4(f,0);
            }
            
            ////////////////////////////

            torqueVector = real(-1.0)*torqueVector;

            //Axis
            prefactor = real(0);
            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];
                
                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com2;
                
                real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);
                prefactor+=mi*dot(ui,ui);
            }
            prefactor=real(1.0)/prefactor;
            
            for(int i=0;i<setSize2;i++){
                int index = id2index[set2id2[i+set*setSize2]];

                real  mi = mass[index];
                real3 ri = make_real3(pos[index])-com2;
                
                real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);
                
                const real3 f =  mi*prefactor*cross(torqueVector,ui);

                force[index] += make_real4(f,0);
            }
        }
    }
    
    class ConstantTorqueBtwCOM: public Interactor{

        protected:

            int nSets;
            
            int setSize1;
            int setSize2;
            
            thrust::device_vector<int> set2id1;
            thrust::device_vector<int> set2id2;
            
            thrust::device_vector<real> T;
            
            Box box;

        public:

            struct Parameters{};
            
            ConstantTorqueBtwCOM(std::shared_ptr<System>       sys,
                                 std::shared_ptr<ParticleData>  pd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 int setSize1,
                                 int setSize2,
                                 int nSets,
                                 thrust::host_vector<int> set2id1,
                                 thrust::host_vector<int> set2id2,
                                 thrust::host_vector<real> T,
                                 Parameters par):Interactor(pd,pg,sys,std::string("ConstantTorqueBtwCOM")),
                                                 setSize1(setSize1),
                                                 setSize2(setSize2),
                                                 nSets(nSets),
                                                 set2id1(set2id1),
                                                 set2id2(set2id2),
                                                 T(T){}
            
            void sum(Computables comp,cudaStream_t st) override {

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                int Nthreads = 128;
                int Nblocks=nSets/Nthreads + ((nSets%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    ConstantTorqueBtwCOM_ns::computeConstantTorqueBtwCOM_TorqueKernel<<<Nblocks, Nthreads, 0, st>>>(pd->getIdOrderedIndices(access::location::gpu),
                                                                                                                    mass,
                                                                                                                    pos,
                                                                                                                    force,
                                                                                                                    setSize1, 
                                                                                                                    setSize2, 
                                                                                                                    nSets, 
                                                                                                                    thrust::raw_pointer_cast(set2id1.data()), 
                                                                                                                    thrust::raw_pointer_cast(set2id2.data()), 
                                                                                                                    thrust::raw_pointer_cast(T.data()),  
                                                                                                                    box);  

                }
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
}}}}


#endif
