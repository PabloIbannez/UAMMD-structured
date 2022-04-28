#ifndef __INTEGRATOR_BASIC__
#define __INTEGRATOR_BASIC__

#include "Integrator/Integrator.cuh"
         
namespace uammd{
namespace structured{

    namespace IntegratorBasic_ns{

        template<class PropItr,class PropType>
        void loadProperty(shared_ptr<ParticleData>  pd,
                          shared_ptr<ParticleGroup> pg,
                          PropItr&             propItr,
                          PropType           propValue){
            
            auto id = pd->getId(access::location::cpu, 
                                access::mode::read);
            
            auto groupIndex  = pg->getIndexIterator(access::location::cpu);
            auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            fori(0,pg->getNumberParticles()){
                int index = sortedIndex[id[groupIndex[i]]];
                propItr[i] = propValue;
            }

        }

        void loadFrictionConstant(shared_ptr<ParticleData>  pd,
                                  shared_ptr<ParticleGroup> pg,
                                  real frictionConstant){

            auto fricConst = pd->getFrictionConstant(access::location::cpu, 
                                                     access::mode::write);

            loadProperty<decltype(fricConst),real>(pd,pg,fricConst,frictionConstant);
        }
        
        void generateVelocity(shared_ptr<ParticleData>  pd,
                              shared_ptr<ParticleGroup> pg,
                              shared_ptr<System>       sys,
                              int  N,
                              real kBT,
                              std::string name,
                              cudaStream_t stream){

            if(pd->isVelAllocated()){
                sys->log<System::WARNING>("[%s] Velocity overwritten!", name.c_str());
            }

            auto groupIterator = pg->getIndexIterator(access::location::gpu);

            auto mass = pd->getMass(access::location::gpu, access::mode::read);

            auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
            auto vel = pd->getVel(access::location::gpu, access::mode::readwrite);

            uint seed_v = sys->rng().next32();

            real* mass_ptr    = mass.raw();

            real4* pos_ptr     = pos.raw();
            real3* vel_ptr     = vel.raw();

            real sigma_kBT = sqrt(kBT);  
            thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                    [=] __host__ __device__ (int index){const real  sigma = sigma_kBT/sqrt(mass_ptr[index]);
                    Saru rng_s(index, seed_v);
                    const real3 noise = make_real3(rng_s.gf(0, sigma), rng_s.gf(0, sigma).x);
                    vel_ptr[index]=noise;});

            //CudaSafeCall(cudaStreamSynchronize(stream));

        }

    }

    class IntegratorBasic: public Integrator{
        
        protected:
            
            uint seed;

            cudaStream_t stream;
            uint steps=0;

            bool constrained= false;
            shared_ptr<Constraint::Constraint> constraint;

        public:

            struct Parameters{

                Box box;

                real dt = 0;

                int  stopTransRotSteps=0;
            };

        protected:
            
            Box box;
                
            real dt;

            int  stopTransRotSteps;

        public:
            
            Parameters inputFileToParam(InputFile& in){

                Parameters param;

                real3 boxSizeBuffer;

                in.getOption("boxSize",InputFile::Required)>>boxSizeBuffer.x
                                                           >>boxSizeBuffer.y
                                                           >>boxSizeBuffer.z;
                
                param.box = Box(boxSizeBuffer);
                
                in.getOption("dt",InputFile::Required)>>param.dt;
                
                in.getOption("stopTransRotSteps",InputFile::Optional)
                                                 >>param.stopTransRotSteps;

                return param;
            }

        IntegratorBasic(shared_ptr<ParticleData>  pd,
                        shared_ptr<ParticleGroup> pg,
                        shared_ptr<System>       sys,
                        Parameters param,
                        std::string name,
                        cudaStream_t stream):Integrator(pd, pg, sys, name),
                                             box(param.box),
                                             dt(param.dt),
                                             stopTransRotSteps(param.stopTransRotSteps),
                                             stream(stream){
                        
                            sys->log<System::MESSAGE>("[%s] Box: (%f,%f,%f)", name.c_str(),
                                                                              box.boxSize.x,
                                                                              box.boxSize.y,
                                                                              box.boxSize.z);
                            
                            sys->log<System::MESSAGE>("[%s] Time step: %f",name.c_str(), 
                                                                           dt);

                            sys->log<System::MESSAGE>("[%s] stopTransRotSteps: %i", name.c_str(),
                                                                                    stopTransRotSteps);

                            sys->rng().next32();
                            sys->rng().next32();
                            seed = sys->rng().next32();
                
                        }

            void setStep(uint nSteps){steps=nSteps;}

            Box  getBox(){return box;}
            real getTimeStep(){return dt;}

            void setConstraint(std::shared_ptr<Constraint::Constraint> newConstraint){
                if(constrained==false){
                    this->constraint = newConstraint;
                    constrained = true;
                } else {
                    sys->log<System::CRITICAL>("[%s] Trying to add constraint %s, but constraint %s has yet been added.", name.c_str(),
                                                                                                                          newConstraint->getName(),
                                                                                                                          constraint->getName());
                }
            }
            
            template<class UNITS>
            void applyUnits(){
                dt = dt*UNITS::TO_INTERNAL_TIME;
                sys->log<System::MESSAGE>("[%s] Time step (after units): %f", this->name.c_str(), dt);
            }
            
            void resetForce(cudaStream_t st){
                
                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), force.begin(), force.end(), make_real4(0));
                
                CudaCheckError();
            }

            void resetForce(){
                resetForce(stream);
            }
            
            void resetEnergy(cudaStream_t st){
                
                auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), real(0));
                
                CudaCheckError();
            }
            
            void resetEnergy(){
                resetEnergy(stream);
            }
            
            void resetVirial(cudaStream_t st){
                
                auto virial = pd->getVirial(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), virial.begin(), virial.end(), tensor3(0));
                
                CudaCheckError();
            }
            
            void resetVirial(){
                resetVirial(stream);
            }
            
            void sumForce() {
                this->sumForce(stream);
            }
            
            void sumForce(cudaStream_t st){
                for(auto forceComp: interactors) forceComp->sum({.force =true,
                                                                 .energy=false,
                                                                 .virial=false,},st);
                CudaCheckError();
            }
            
            void sumVirial() {
                this->sumVirial(stream);
            }
            
            void sumVirial(cudaStream_t st){
                for(auto virialComp: interactors) virialComp->sum({.force =false,
                                                                   .energy=false,
                                                                   .virial=true,},st);
                CudaCheckError();
            }
            
            void sumEnergy() override {
                this->sumEnergy(stream);
            }

            void sumEnergy(cudaStream_t st){
                for(auto energyComp: interactors) energyComp->sum({.force =false,
                                                                   .energy=true,
                                                                   .virial=false,},st);
                CudaCheckError();
            }
            
            void update(){
                for(auto upd: updatables){
                    upd->updateSimulationTime(steps*dt);
                    upd->updateTimeStep(dt);
                    upd->updateBox(box);
                }
            }
            

    };

    class IntegratorBasicNVT: public IntegratorBasic {
        
        public:

            struct Parameters: public IntegratorBasic::Parameters{
                real kB = 1;
                real T  = 0;
            };

        protected:
            
            real kB = 1;
            real T  = 0;

        public:
            
            Parameters inputFileToParam(InputFile& in){
                
                Parameters param;
                static_cast<IntegratorBasic::Parameters&>(param) = IntegratorBasic::inputFileToParam(in); 
            
                in.getOption("kBoltzmann",InputFile::Optional)
                              >>param.kB;
                
                in.getOption("T",InputFile::Required)
                              >>param.T;
                
                return param;
            }

        IntegratorBasicNVT(shared_ptr<ParticleData>  pd,
                           shared_ptr<ParticleGroup> pg,
                           shared_ptr<System>       sys,
                           Parameters param,
                           std::string name,
                           cudaStream_t stream):IntegratorBasic(pd, pg, sys,param,name,stream),
                                                kB(param.kB), 
                                                T(param.T){
                        
                sys->log<System::MESSAGE>("[%s] kB: %f",name.c_str(),kB);
                sys->log<System::MESSAGE>("[%s] T: %f" ,name.c_str(),T);
                                                
            }
            
            template<class UNITS>
            void applyUnits(){
                IntegratorBasic::applyUnits<UNITS>();
                
                kB = UNITS::KBOLTZ;
                sys->log<System::MESSAGE>("[%s] kB (after units): %f", this->name.c_str(),kB);
                sys->log<System::MESSAGE>("[%s] kBT (after units): %f", this->name.c_str(),kB*T);
            }
            
            void update(){
                IntegratorBasic::update();
                for(auto upd: updatables){
                    upd->updateTemperature(T);
                }
            }
    };

}}

#endif

