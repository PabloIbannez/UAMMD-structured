#include "Integrator/Hydro/ICM.cuh"

#ifndef __ICM__
#define __ICM__

namespace uammd{
namespace structured{

    namespace Hydrodynamics{

    class ICM: public uammd::Hydro::ICM{

            real kB = 1;
            real  T = 0;

        public:

            struct Parameters: public uammd::Hydro::ICM::Parameters{
                real kB = 1;
                real  T = 0;
            };

        private:

            Parameters inputFileToParam(InputFile& in){
                
                Parameters param;
                
                real3 boxSizeBuffer;
                
                in.getOption("boxSize",InputFile::Required)>>boxSizeBuffer.x
                                                           >>boxSizeBuffer.y
                                                           >>boxSizeBuffer.z;
                
                param.box = Box(boxSizeBuffer);
                
                in.getOption("dt",InputFile::Required)
                              >>param.dt;
                
                in.getOption("kBoltzmann",InputFile::Optional)
                              >>param.kB;
                in.getOption("T",InputFile::Required)
                              >>param.T;

                param.temperature = param.kB*param.T; 
                
                in.getOption("viscosity",InputFile::Required)
                              >>param.viscosity;
                in.getOption("density",InputFile::Required)
                              >>param.density;
                
                in.getOption("hydrodynamicRadius",InputFile::Required)
                              >>param.hydrodynamicRadius;


                return param;
            }

        public:

            ICM(shared_ptr<ParticleData>  pd,
                shared_ptr<ParticleGroup> pg,
                shared_ptr<System>       sys,
                InputFile in                ,
                cudaStream_t stream):uammd::Hydro::ICM(pd,sys,inputFileToParam(in)){

                    Parameters param = inputFileToParam(in);

                    T  = param.T;
                    kB = param.kB;

                    st=stream;
                            
                    sys->log<System::MESSAGE>("[%s] Box: (%f,%f,%f)", name.c_str(),
                                                                      box.boxSize.x,
                                                                      box.boxSize.y,
                                                                      box.boxSize.z);
                    
                    sys->log<System::MESSAGE>("[%s] Time step: %f",name.c_str(), 
                                                                   dt);
                
                    sys->log<System::MESSAGE>("[%s] kB: %f",name.c_str(),kB);
                    sys->log<System::MESSAGE>("[%s] T: %f" ,name.c_str(),T);
                    
            }
            
            void resetForce(cudaStream_t st){
                
                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), force.begin(), force.end(), make_real4(0));
                
                CudaCheckError();
            }

            void resetForce(){
                resetForce(st);
            }
            
            void resetEnergy(cudaStream_t st){
                
                auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), real(0));
                
                CudaCheckError();
            }
            
            void resetEnergy(){
                resetEnergy(st);
            }
            
            void resetVirial(cudaStream_t st){
                
                auto virial = pd->getVirial(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(st), virial.begin(), virial.end(), tensor3(0));
                
                CudaCheckError();
            }
            
            void resetVirial(){
                resetVirial(st);
            }
            
            void sumForce() {
                this->sumForce(st);
            }
            
            void sumForce(cudaStream_t st){
                for(auto forceComp: interactors) forceComp->sum({.force =true,
                                                                 .energy=false,
                                                                 .virial=false,},st);
                CudaCheckError();
            }
            
            void sumVirial() {
                this->sumVirial(st);
            }
            
            void sumVirial(cudaStream_t st){
                for(auto virialComp: interactors) virialComp->sum({.force =false,
                                                                   .energy=false,
                                                                   .virial=true,},st);
                CudaCheckError();
            }
            
            void sumEnergy() override {
                this->sumEnergy(st);
            }

            void sumEnergy(cudaStream_t st){
                for(auto energyComp: interactors) energyComp->sum({.force =false,
                                                                   .energy=true,
                                                                   .virial=false,},st);
                CudaCheckError();
            }
            
            template<class UNITS>
            void applyUnits(){
                
                dt = dt*UNITS::TO_INTERNAL_TIME;
                sys->log<System::MESSAGE>("[%s] Time step (after units): %f", this->name.c_str(), dt);
                
                kB = UNITS::KBOLTZ;
                sys->log<System::MESSAGE>("[%s] kB (after units): %f", this->name.c_str(),kB);

                temperature = T*kB;
                sys->log<System::MESSAGE>("[%s] kBT (after units): %f", this->name.c_str(),temperature);
                
                viscosity = viscosity/UNITS::TO_INTERNAL_TIME;
                sys->log<System::MESSAGE>("[%s] Viscosity (after units): %f", this->name.c_str(), viscosity);
            }
            
            void update(){
                for(auto upd: updatables){
                    upd->updateSimulationTime(step*dt);
                    upd->updateTimeStep(dt);
                    upd->updateBox(box);
                    upd->updateTemperature(temperature);
                    upd->updateViscosity(viscosity);
                }
            }



    };


    }
}}

#endif
