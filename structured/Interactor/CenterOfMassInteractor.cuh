#ifndef __COM_INTERACTOR__
#define __COM_INTERACTOR__

namespace uammd{
namespace structured{
namespace Interactor{
    
    namespace ExternalCOM_ns{

        void __global__ computeExternalCOMForceKernel(real*  mass,
                                                      real4* pos,
                                                      real4* force,
                                                      int numberParticles1,
			                                          ParticleGroup::IndexIterator group1Iterator,
                                                      Box box,
                                                      real M1,
                                                      real3 F){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= numberParticles1) return;

            int index = group1Iterator[i];
            const real3 f =  (mass[index]/M1)*F;
            force[index] += make_real4(f,0);

        }
    }

    class ExternalCOM: public Interactor{

        private:

            int N1;
            
            std::shared_ptr<ParticleGroup> group1;
            
            real M1;
            
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;

            Box box;

            //State 

            real3 F;

        public:

            struct Parameters{};

            ExternalCOM(std::shared_ptr<ParticleGroup> pg,
                        std::shared_ptr<ParticleGroup> group1,
                        Parameters par):Interactor(pg,std::string("ExternalCOM")),
                                        group1(group1)
            {
                N1=group1->getNumberParticles();

                M1 = Measures::totalMass(group1,0);

                this->setState({0.0,0.0,0.0});

                cudaDeviceSynchronize();
            }

            ~ExternalCOM(){
                cudaFree(cubReductionTempStorage);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                    
                auto g1Iter = group1->getIndexIterator(access::location::gpu);

                real*  mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos  = pd->getPos(access::location::gpu, access::mode::read).raw();

                int blocksize = 256;
                int Nthreads = blocksize<N1?blocksize:N1;
                int Nblocks=N1/Nthreads + ((N1%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    ExternalCOM_ns::computeExternalCOMForceKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                pos,
                                                                                                force,
                                                                                                N1,
			                                                                                    g1Iter,
                                                                                                box,
                                                                                                M1,
                                                                                                F);  

                }
            }
            
            void  setState(real3 newForce){F=newForce;}
            real3 getState(){return F;}
            
            void measure(std::ofstream& outputFile){
                outputFile << this->getState() << std::endl;
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
            
    };
    
    namespace HarmonicCOM_ns{

        void __global__ computeHarmonicCOMForceKernel(real*  mass,
                                                      real4* pos,
                                                      real4* force,
                                                      int numberParticles1,
                                                      int numberParticles2,
			                                          ParticleGroup::IndexIterator group1Iterator,
			                                          ParticleGroup::IndexIterator group2Iterator,
                                                      Box box,
                                                      real M1,
                                                      real M2,
                                                      real4* com1M1,
                                                      real4* com2M2,
                                                      real K,
                                                      real r0){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= (numberParticles1+numberParticles2)) return;

            const real3 r12 = box.apply_pbc(make_real3(*(com2M2))/M2-make_real3(*(com1M1))/M1);
            
            const real  r2  = dot(r12,r12);
            const real invr = rsqrt(r2);
            
            //Force over 1
            real3 F12 =  Potentials::CommonPotentials::Harmonic::force(r12,r2,K,r0);
            //real3 F12 =  K*(real(1.0)-r0*invr)*r12;

            if(r0 == real(0.0) and r2 == real(0.0)){
                real3 F12 = make_real3(0.0);
            }

            if(i < numberParticles1){
                int index = group1Iterator[i];
                const real3 f =  (mass[index]/M1)*F12;
                force[index] += make_real4(f,0);
            } else {
                i = i - numberParticles1;
                int index = group2Iterator[i];
                const real3 f = -(mass[index]/M2)*F12;
                force[index] += make_real4(f,0);
            }

        }
        
        void __global__ computeHarmonicCOMEnergyKernel(real*  mass,
                                                       real4* pos,
                                                       real*  energy,
                                                       int numberParticles1,
                                                       int numberParticles2,
			                                           ParticleGroup::IndexIterator group1Iterator,
			                                           ParticleGroup::IndexIterator group2Iterator,
                                                       Box box,
                                                       real M1,
                                                       real M2,
                                                       real4* com1M1,
                                                       real4* com2M2,
                                                       real K,
                                                       real r0){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= (numberParticles1+numberParticles2)) return;

            const real3 r12 = box.apply_pbc(make_real3(*(com2M2))/M2-make_real3(*(com1M1))/M1);
            
            const real r2 = dot(r12,r12);
            //const real r = sqrt(dot(r12,r12));
            
            //Energy
            real e =  Potentials::CommonPotentials::Harmonic::energy(r12,r2,K,r0);
            //real e =  real(0.5)*K*(r-r0)*(r-r0);
            
            if(i < numberParticles1){
                int index = group1Iterator[i];
                energy[index]+=real(0.5)*(mass[index]/M1)*e;
            } else {
                i = i - numberParticles1;
                int index = group2Iterator[i];
                energy[index]+=real(0.5)*(mass[index]/M2)*e;
            }

        }
    }

    class HarmonicCOM: public Interactor{

        private:

            int N1;
            int N2;
            
            std::shared_ptr<ParticleGroup> group1;
            std::shared_ptr<ParticleGroup> group2;
            
            real M1;
            real M2;

            thrust::device_vector<real4> com1M1;
            thrust::device_vector<real4> com2M2;
            
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;

            Box box;

            //State 

            real r0;

        public:

            struct Parameters{
                real K;
            };

        private:

            real K;

        public:
            
            HarmonicCOM(std::shared_ptr<ParticleGroup> pg,
                        std::shared_ptr<ParticleGroup> group1,
                        std::shared_ptr<ParticleGroup> group2,
                        Parameters par):Interactor(pg,std::string("HarmonicCOM")),
                                        group1(group1),
                                        group2(group2),
                                        K(par.K)
            {
                N1=group1->getNumberParticles();
                N2=group2->getNumberParticles();

                M1 = Measures::totalMass(group1,0);
                M2 = Measures::totalMass(group2,0);

                com1M1.resize(1);
                com2M2.resize(1);

                cubReductionTempStorage     = NULL;
                cubReductionTempStorageSize = 0;
        
                {
                    real*  mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                    real4* pos  = pd->getPos(access::location::gpu, access::mode::read).raw();
                    
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);
        
                    auto g1Iter = group1->getIndexIterator(access::location::gpu);
                    auto g2Iter = group2->getIndexIterator(access::location::gpu);
                    
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g1Iter, mWs), //It does not matter if 1 or 2 at init 
                                           thrust::raw_pointer_cast(com1M1.data()), 
                                           std::max(N1,N2));

                    cudaMalloc(&cubReductionTempStorage, cubReductionTempStorageSize);
                        
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g1Iter, mWs), 
                                           thrust::raw_pointer_cast(com1M1.data()), 
                                           N1);
                    
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g2Iter, mWs), 
                                           thrust::raw_pointer_cast(com2M2.data()), 
                                           N2);
                }
                
                this->setState(this->getDistance());

                cudaDeviceSynchronize();
            }

            ~HarmonicCOM(){
                cudaFree(cubReductionTempStorage);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                    
                auto g1Iter = group1->getIndexIterator(access::location::gpu);
                auto g2Iter = group2->getIndexIterator(access::location::gpu);

                real*  mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos  = pd->getPos(access::location::gpu, access::mode::read).raw();

                {
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);

                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                            cubReductionTempStorageSize,
                            thrust::make_transform_iterator(g1Iter, mWs), 
                            thrust::raw_pointer_cast(com1M1.data()), 
                            N1,st);

                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                            cubReductionTempStorageSize,
                            thrust::make_transform_iterator(g2Iter, mWs), 
                            thrust::raw_pointer_cast(com2M2.data()), 
                            N2,st);
                }

                int N = N1 + N2;

                int blocksize = 256;
                int Nthreads = blocksize<N?blocksize:N;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    HarmonicCOM_ns::computeHarmonicCOMForceKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                pos,
                                                                                                force,
                                                                                                N1,N2,
			                                                                                    g1Iter,g2Iter,
                                                                                                box,
                                                                                                M1,M2,
                                                                                                thrust::raw_pointer_cast(com1M1.data()), 
                                                                                                thrust::raw_pointer_cast(com2M2.data()), 
                                                                                                K,
                                                                                                r0);  

                }
                
                if(comp.energy == true){
                    
                    auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

                    HarmonicCOM_ns::computeHarmonicCOMEnergyKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                 pos,
                                                                                                 energy,
                                                                                                 N1,N2,
			                                                                                     g1Iter,g2Iter,
                                                                                                 box,
                                                                                                 M1,M2,
                                                                                                 thrust::raw_pointer_cast(com1M1.data()), 
                                                                                                 thrust::raw_pointer_cast(com2M2.data()), 
                                                                                                 K,
                                                                                                 r0);  




                }

                //CudaSafeCall(cudaStreamSynchronize(st));
            }
            
            void   setState(double newDistance){r0=newDistance;}
            double getState(){return r0;}
            
            real3 getCOM1(){
                return make_real3(com1M1[0]/M1);
            }
            
            real3 getCOM2(){
                return make_real3(com2M2[0]/M2);
            }

            real getDistance(){
                real3 r12 = box.apply_pbc(this->getCOM2()-this->getCOM1());
                return sqrt(dot(r12,r12));
            }
            
            void measure(std::ofstream& outputFile){
                real dr = this->getDistance()-this->getState();
                outputFile << this->getState()     << " " 
                           << this->getDistance()  << " " 
                           << real(0.5)*K*dr*dr << " " << K*dr << std::endl;
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
            
    };
    
    namespace HarmonicFixedCOM_ns{

        void __global__ computeHarmonicFixedCOMForceKernel(real* mass,
                                                           real4* pos,
                                                           real4* force,
                                                           int numberParticles,
			                                               ParticleGroup::IndexIterator groupIterator,
                                                           Box box,
                                                           real M,
                                                           real4* comM,
                                                           real3 K,
                                                           real3 fixedPoint){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= numberParticles) return;

            const real3 r1f = box.apply_pbc(fixedPoint-make_real3(*(comM))/M);

            //Force over 1
            real3 F1f =  Potentials::CommonPotentials::HarmonicAnisotropic::force(r1f,
                                                                                  K,
                                                                                  make_real3(0.0));
            int index = groupIterator[i];
            const real3 f = (mass[index]/M)*F1f;
            
            force[index] += make_real4(f,0);
        }
    }
    
    //Fixed to a point
    class HarmonicFixedCOM: public Interactor{

        private:

            int N;
            
            real M;

            thrust::device_vector<real4> comM;
            
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;
                
            Box box;
                                               
            //State

            real3 fixedPoint;

        public:

            struct Parameters{
                real3 K;
            };

        private:

            real3 K;

        public:
            
            HarmonicFixedCOM(std::shared_ptr<ParticleGroup> pg,
                             Parameters par):Interactor(pg,std::string("HarmonicFixedCOM")),
                                             K(par.K)
            {
                N = pg->getNumberParticles();

                M = Measures::totalMass(pg,0);

                comM.resize(1);

                cubReductionTempStorage     = NULL;
                cubReductionTempStorageSize = 0;
        
                {
                    real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                    real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
                    
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);
        
                    auto gIter = pg->getIndexIterator(access::location::gpu);
                    
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(gIter, mWs), 
                                           thrust::raw_pointer_cast(comM.data()), 
                                           N);

                    cudaMalloc(&cubReductionTempStorage, cubReductionTempStorageSize);
                        
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(gIter, mWs), 
                                           thrust::raw_pointer_cast(comM.data()), 
                                           N);

                    this->setState(this->getCOM());
                }

                cudaDeviceSynchronize();
            }
            
            ~HarmonicFixedCOM(){
                cudaFree(cubReductionTempStorage);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                    
                auto gIter = pg->getIndexIterator(access::location::gpu);
                    
                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
                
                {
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);
        
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(gIter, mWs), 
                                           thrust::raw_pointer_cast(comM.data()), 
                                           N,st);
                }

                int blocksize = 256;
                int Nthreads  = blocksize<N?blocksize:N;
                int Nblocks   = N/Nthreads + ((N%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
                    
                    HarmonicFixedCOM_ns::computeHarmonicFixedCOMForceKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                          pos,
                                                                                                          force,
                                                                                                          N,
			                                                                                              gIter,
                                                                                                          box,
                                                                                                          M,
                                                                                                          thrust::raw_pointer_cast(comM.data()), 
                                                                                                          K,
                                                                                                          fixedPoint);  

                }
                //CudaSafeCall(cudaStreamSynchronize(st));
            }
            
            real3 getCOM(){
                return make_real3(comM[0]/M);
            }
            real getDistance(){
                real3 r12 = box.apply_pbc(fixedPoint-this->getCOM());
                return sqrt(dot(r12,r12));
            }
            
            void measure(std::ofstream& outputFile){
                real3 r12 = box.apply_pbc(fixedPoint-this->getCOM());
                real e  = real(0.5)*K.x*(r12.x-fixedPoint.x)*(r12.x-fixedPoint.x)+
                          real(0.5)*K.y*(r12.y-fixedPoint.y)*(r12.y-fixedPoint.y)+
                          real(0.5)*K.z*(r12.z-fixedPoint.z)*(r12.z-fixedPoint.z);
                real3 f = {K.x*(r12.x-fixedPoint.x),
                           K.y*(r12.y-fixedPoint.y),
                           K.z*(r12.z-fixedPoint.z)};
                outputFile << this->getDistance()  << " " 
                           << e << " " << f << std::endl;
            }

            void  setState(real3 newfixedPoint){fixedPoint=newfixedPoint;}
            real3 getState(){return fixedPoint;}
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
    
    namespace ConstantForceCOM_ns{

        void __global__ computeConstantForceCOMForceKernel(real* mass,
                                                           real4* pos,
                                                           real4* force,
                                                           int numberParticles1,
                                                           int numberParticles2,
			                                               ParticleGroup::IndexIterator group1Iterator,
			                                               ParticleGroup::IndexIterator group2Iterator,
                                                           Box box,
                                                           real M1,
                                                           real M2,
                                                           real4* com1M1,
                                                           real4* com2M2,
                                                           real F){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            if(i >= (numberParticles1+numberParticles2)) return;

            const real3 r12 = box.apply_pbc(make_real3(*(com2M2))/M2-make_real3(*(com1M1))/M1);
            
            const real  r2  = dot(r12,r12);
            const real invr = rsqrt(r2);
            
            //Force over 1
            real3 F12 =  -F*r12*invr;

            if(r2 == real(0.0)){
                real3 F12 = make_real3(0.0);
            }

            if(i < numberParticles1){
                int index = group1Iterator[i];
                const real3 f =  (mass[index]/M1)*F12;
                force[index] += make_real4(f,0);
            } else {
                i = i - numberParticles1;
                int index = group2Iterator[i];
                const real3 f = -(mass[index]/M2)*F12;
                force[index] += make_real4(f,0);
            }

        }
    }

    class ConstantForceCOM: public Interactor{

        private:

            int N1;
            int N2;
            
            std::shared_ptr<ParticleGroup> group1;
            std::shared_ptr<ParticleGroup> group2;
                
            real M1;
            real M2;

            thrust::device_vector<real4> com1M1;
            thrust::device_vector<real4> com2M2;
            
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;

            Box box;
            
            real F;

        public:

            struct Parameters{};
            
            ConstantForceCOM(std::shared_ptr<ParticleGroup> pg,
                             std::shared_ptr<ParticleGroup> group1,
                             std::shared_ptr<ParticleGroup> group2,
                             Parameters par):Interactor(pg,std::string("ConstantForceCOM")),
                                             group1(group1),
                                             group2(group2),
                                             F(0)
            {
                N1=group1->getNumberParticles();
                N2=group2->getNumberParticles();

                M1 = Measures::totalMass(group1,0);
                M2 = Measures::totalMass(group2,0);

                com1M1.resize(1);
                com2M2.resize(1);

                cubReductionTempStorage     = NULL;
                cubReductionTempStorageSize = 0;
        
                {
                    real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                    real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
                    
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);
        
                    auto g1Iter = group1->getIndexIterator(access::location::gpu);
                    auto g2Iter = group2->getIndexIterator(access::location::gpu);
                    
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g1Iter, mWs), //It does not matter if 1 or 2 at init 
                                           thrust::raw_pointer_cast(com1M1.data()), 
                                           std::max(N1,N2));

                    cudaMalloc(&cubReductionTempStorage, cubReductionTempStorageSize);
                        
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g1Iter, mWs), 
                                           thrust::raw_pointer_cast(com1M1.data()), 
                                           N1);
                    
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::make_transform_iterator(g2Iter, mWs), 
                                           thrust::raw_pointer_cast(com2M2.data()), 
                                           N2);
                }

                cudaDeviceSynchronize();
            }
            
            ~ConstantForceCOM(){
                cudaFree(cubReductionTempStorage);
            }
            
            void sum(Computables comp,cudaStream_t st) override {

                auto g1Iter = group1->getIndexIterator(access::location::gpu);
                auto g2Iter = group2->getIndexIterator(access::location::gpu);

                real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                {
                    Measures::MeasuresTransforms::weightedSum<real4> mWs(mass,pos);

                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                            cubReductionTempStorageSize,
                            thrust::make_transform_iterator(g1Iter, mWs), 
                            thrust::raw_pointer_cast(com1M1.data()), 
                            N1,st);

                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                            cubReductionTempStorageSize,
                            thrust::make_transform_iterator(g2Iter, mWs), 
                            thrust::raw_pointer_cast(com2M2.data()), 
                            N2,st);
                }

                int N = N1 + N2;

                int blocksize = 256;
                int Nthreads = blocksize<N?blocksize:N;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                if(comp.force == true){

                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    ConstantForceCOM_ns::computeConstantForceCOMForceKernel<<<Nblocks, Nthreads, 0, st>>>(mass,
                                                                                                          pos,
                                                                                                          force,
                                                                                                          N1,N2,
			                                                                                              g1Iter,g2Iter,
                                                                                                          box,
                                                                                                          M1,M2,
                                                                                                          thrust::raw_pointer_cast(com1M1.data()), 
                                                                                                          thrust::raw_pointer_cast(com2M2.data()), 
                                                                                                          F);  

                }
            }
            
            void   setState(double newConstantForce){F=newConstantForce;}
            double getState(){return F;}
            
            real3 getCOM1(){
                return make_real3(com1M1[0]/M1);
            }
            
            real3 getCOM2(){
                return make_real3(com2M2[0]/M2);
            }

            real getDistance(){
                real3 r12 = box.apply_pbc(this->getCOM2()-this->getCOM1());
                return sqrt(dot(r12,r12));
            }
            
            void measure(std::ofstream& outputFile){
                real r = this->getDistance();
                outputFile << this->getDistance()  << " " 
                           << -F*r << " " << -F << std::endl;
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
            }
    };
}}}


#endif
