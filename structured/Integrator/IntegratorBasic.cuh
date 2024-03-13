#ifndef __INTEGRATOR_BASIC__
#define __INTEGRATOR_BASIC__

#include "Integrator/Integrator.cuh"
#include "IntegratorUtils.cuh"

namespace uammd{
namespace structured{

    class IntegratorBasic: public Integrator{

        protected:

            cudaStream_t stream;

            std::shared_ptr<GlobalData>           gd;

            real dt;

        public:

            IntegratorBasic(std::shared_ptr<GlobalData>           gd,
                            std::shared_ptr<ParticleGroup>        pg,
                            DataEntry& data,
                            std::string name):Integrator(pg,name),gd(gd){

                System::log<System::DEBUG>("[IntegratorBasic] Created integrator \"%s\"",name.c_str());

                stream = gd->getSystem()->getCudaStream();

                if(data.isParameterAdded("timeStep")){
                    dt = data.getParameter<real>("timeStep");
                    gd->getFundamental()->setTimeStep(dt);
                    System::log<System::DEBUG>("[IntegratorBasic] (%s) Reading time step (timeStep=%f) from input file",name.c_str(),dt);
                }else{
                    dt = gd->getFundamental()->getTimeStep();
                }

                System::log<System::MESSAGE>("[IntegratorBasic] (%s) Time step: %f",name.c_str(),dt);
                System::log<System::MESSAGE>("[IntegratorBasic] (%s) Current step: %llu",name.c_str(),gd->getFundamental()->getCurrentStep());

            }

            void resetEnergy(){

                auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite);
                thrust::fill(thrust::cuda::par.on(stream), energy.begin(), energy.end(), real(0));

                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void resetForce(){

                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);
                thrust::fill(thrust::cuda::par.on(stream), force.begin(), force.end(), make_real4(0));

                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

           void resetTorque(){

                auto torque = pd->getTorque(access::location::gpu, access::mode::readwrite);
                thrust::fill(thrust::cuda::par.on(stream), torque.begin(), torque.end(), make_real4(0));

                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void updateEnergy(){
                for(auto energyComp: interactors) energyComp->sum({.force =false,
                                                                   .energy=true,
                                                                   .virial=false,
                                                                   .stress=false},stream);
                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void updateForce(bool computeMagneticField = false) {
                for(auto forceComp: interactors) forceComp->sum({.force =true,
                                                                 .energy=false,
                                                                 .virial=false,
                                                                 .stress=false,
		                                                 .magneticField=computeMagneticField},
		                                                 stream);
                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            virtual void forwardTime() override = 0;

    };

    class IntegratorBasicNVT: public IntegratorBasic {

        protected:

            real temperature;
            real kBT;

        public:

            IntegratorBasicNVT(std::shared_ptr<GlobalData>           gd,
                               std::shared_ptr<ParticleGroup>        pg,
                               DataEntry& data,
                               std::string name):IntegratorBasic(gd,pg,data,name){

                System::log<System::DEBUG>("[IntegratorBasicNVT] Created integrator \"%s\"",name.c_str());

                //Look for the temperature in the input file
                if(data.isParameterAdded("temperature")){
                    temperature = data.getParameter<real>("temperature");
                    System::log<System::DEBUG>("[IntegratorBasicNVT] (%s) Reading temperature (T=%f) from input file",name.c_str(),temperature);
                }else{
                    temperature = gd->getEnsemble()->getTemperature();
                }
                kBT         = gd->getUnits()->getBoltzmannConstant()*temperature;

                System::log<System::MESSAGE>("[IntegratorBasicNVT] (%s) Temperature: %f",name.c_str(),temperature);
                System::log<System::MESSAGE>("[IntegratorBasicNVT] (%s) kBT (kB %f): %f",name.c_str(),gd->getUnits()->getBoltzmannConstant(),kBT);
            }

            virtual void forwardTime() override = 0;

    };

}}


#endif

