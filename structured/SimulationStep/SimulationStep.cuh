#ifndef __SIMULATION_STEP__
#define __SIMULATION_STEP__

namespace uammd{
namespace structured{
namespace SimulationStep{

class SimulationStepBase{

    protected:

        std::shared_ptr<ParticleGroup>              pg;
        std::shared_ptr<IntegratorManager>  integrator;
        std::shared_ptr<ForceField>                 ff;

        std::shared_ptr<Topology>       topology;

        std::shared_ptr<ExtendedSystem>      sys;
        std::shared_ptr<GlobalData>           gd;
        std::shared_ptr<ExtendedParticleData> pd;

        std::string name;

        bool initialized = false;

        ullint startStep;
        ullint endStep;
        ullint intervalStep;

        ullint lastStepApplied;

        virtual void init(cudaStream_t st) = 0;
        virtual void applyStep(ullint step, cudaStream_t st) = 0;

        virtual void update(ullint step, cudaStream_t st) {}; // Do nothing by default

        bool isStepApplied(bool force=false){
            int step = this->gd->getFundamental()->getCurrentStep();
            ullint localStep = step - this->startStep;
            if(initialized){
                if(this->intervalStep==0 and !force){return false;}
                if(step>=this->startStep and step<=this->endStep and (step%localStep==0 or force)){
                    return true;
                }
            } else {
                System::log<System::CRITICAL>("[SimulationStep] (%s) SimulationStep not initialized", name.c_str());
            }
        }

    public:

        SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           std::string name):pg(pg),
                                             integrator(integrator),
                                             ff(ff),
                                             name(name){

            this->topology = ff->getTopology();

            sys = topology->getSystem();
            gd  = topology->getGlobalData();
            pd  = topology->getParticleData();

            lastStepApplied = std::numeric_limits<ullint>::max();

        }

        //Construct from struct
        struct Parameters{
            ullint startStep    = 0;
            ullint endStep      = std::numeric_limits<ullint>::max();
            ullint intervalStep = 0;
        };

        SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           Parameters par,
                           std::string name):SimulationStepBase(pg,integrator,ff,name){

            //Read parameters

            startStep    = par.startStep;
            endStep      = par.endStep;
            intervalStep = par.intervalStep;

            //Print read parameters
            System::log<System::MESSAGE>("[SimulationStep] (%s) startStep: %llu", name.c_str(), startStep);
            System::log<System::MESSAGE>("[SimulationStep] (%s) endStep: %llu", name.c_str(), endStep);
            System::log<System::MESSAGE>("[SimulationStep] (%s) intervalStep: %llu", name.c_str(), intervalStep);
        }

        //Construct from DataEntry
        SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           DataEntry& data,
                           std::string name):SimulationStepBase(pg,integrator,ff,name){

            //Read parameters

            startStep    = data.getParameter<ullint>("startStep", 0);
            endStep      = data.getParameter<ullint>("endStep", std::numeric_limits<ullint>::max());
            intervalStep = data.getParameter<ullint>("intervalStep");

            //Print read parameters
            if (startStep != 0){
                System::log<System::MESSAGE>("[SimulationStep] (%s) startStep: %llu", name.c_str(), startStep);
            }
            if (endStep != std::numeric_limits<ullint>::max()){
                System::log<System::MESSAGE>("[SimulationStep] (%s) endStep: %llu", name.c_str(), endStep);
            }
            if (intervalStep != 0){
                System::log<System::MESSAGE>("[SimulationStep] (%s) intervalStep: %llu", name.c_str(), intervalStep);
            }
        }

        ~SimulationStepBase(){}

        std::string getName(){return name;}

        ullint getStartStep()   {return startStep;}
        ullint getEndStep()     {return endStep;}
        ullint getIntervalStep(){return intervalStep;}

        ullint getLastStepApplied(){return lastStepApplied;}

        void tryInit( cudaStream_t st){
            if(!initialized){
                this->init(st);
                initialized=true;
            }
        };

        void tryInit(){
            cudaDeviceSynchronize();
            tryInit(0);
            cudaDeviceSynchronize();
        }


        void tryApplyStep(cudaStream_t st,bool force=false){
            int step = this->gd->getFundamental()->getCurrentStep();
            if(initialized){
                if(this->intervalStep==0 and !force){return;}
                if(step>=this->startStep and step<=this->endStep){
                    this->update(step,st); //This method is called every step
                    if (step%this->intervalStep==0 or force){
                        this->applyStep(step,st);
                        lastStepApplied = step;
                    }
                }
            } else {
                System::log<System::CRITICAL>("[SimulationStep] (%s) SimulationStep not initialized", name.c_str());
            }
        }

        void tryApplyStep(bool force=false){
            cudaDeviceSynchronize();
            tryApplyStep(0,force);
            cudaDeviceSynchronize();
        }

};

class SimulationStepBase_EnergyForceTorque: public SimulationStepBase{

    protected:

      uninitialized_cached_vector<real>  energyTmp;
			uninitialized_cached_vector<real4> forceTmp;
			uninitialized_cached_vector<real4> torqueTmp;

			void copyToTmp(cudaStream_t st){

        System::log<System::DEBUG>("[SimulationStep_EnergyForceTorque] (%s) Copying to tmp", name.c_str());

				int N = this->pd->getNumParticles();

				energyTmp.resize(N);
				forceTmp.resize(N);
				torqueTmp.resize(N);

        {
            auto energy = this->pd->getEnergy(access::location::gpu, access::mode::read);
            thrust::copy(thrust::cuda::par.on(st),
                         energy.begin(),
                         energy.end(),
                         energyTmp.begin());
        }

				{
						auto force = this->pd->getForce(access::location::gpu, access::mode::read);
						thrust::copy(thrust::cuda::par.on(st),
												 force.begin(),
												 force.end(),
												 forceTmp.begin());
				}

				auto torque = this->pd->getTorque(access::location::gpu, access::mode::read);
				thrust::copy(thrust::cuda::par.on(st),
										 torque.begin(),
										 torque.end(),
										 torqueTmp.begin());

				cudaDeviceSynchronize();

			}

			void copyFromTmp(cudaStream_t st){

        System::log<System::DEBUG>("[SimulationStep_EnergyForceTorque] (%s) Copying from tmp", name.c_str());

				int N = this->pd->getNumParticles();

				{
						auto energy = this->pd->getEnergy(access::location::gpu, access::mode::write);
						thrust::copy(thrust::cuda::par.on(st),
												 energyTmp.begin(),
												 energyTmp.end(),
												 energy.begin());
				}

				{
						auto force = this->pd->getForce(access::location::gpu, access::mode::write);
						thrust::copy(thrust::cuda::par.on(st),
												 forceTmp.begin(),
												 forceTmp.end(),
												 force.begin());
				}

				auto torque = this->pd->getTorque(access::location::gpu, access::mode::write);
				thrust::copy(thrust::cuda::par.on(st),
										 torqueTmp.begin(),
										 torqueTmp.end(),
										 torque.begin());

				cudaDeviceSynchronize();

			}

			void setZero(cudaStream_t st){

        System::log<System::DEBUG>("[SimulationStep_EnergyForceTorque] (%s) Setting to zero", name.c_str());

				{
					auto energy = this->pd->getEnergy(access::location::gpu, access::mode::write);
					thrust::fill(thrust::cuda::par.on(st),
							         energy.begin(),
							         energy.end(),
							         real(0));
				}

				{
					auto force = this->pd->getForce(access::location::gpu, access::mode::write);
					thrust::fill(thrust::cuda::par.on(st),
							         force.begin(),
							         force.end(),
							         make_real4(0.0));
				}

				{
					auto torque = this->pd->getTorque(access::location::gpu, access::mode::write);
					thrust::fill(thrust::cuda::par.on(st),
							         torque.begin(),
							         torque.end(),
							         make_real4(0.0));
				}

				cudaDeviceSynchronize();

			}

    public:

        SimulationStepBase_EnergyForceTorque(std::shared_ptr<ParticleGroup>              pg,
                                             std::shared_ptr<IntegratorManager>  integrator,
                                             std::shared_ptr<ForceField>                 ff,
                                             DataEntry& data,
                                             std::string name):SimulationStepBase(pg,integrator,ff,data,name){

        }

        ~SimulationStepBase_EnergyForceTorque(){}

        void tryApplyStep(cudaStream_t st,bool force=false){
            if(this->isStepApplied(force)){
                this->copyToTmp(st);
                SimulationStepBase::tryApplyStep(st,true);
                this->copyFromTmp(st);
            }
        }

        void tryApplyStep(bool force=false){
            if(this->isStepApplied(force)){
                this->copyToTmp(0);
                SimulationStepBase::tryApplyStep(true);
                this->copyFromTmp(0);
            }
        }
};

class SimulationStepBase_EnergyForceTorqueHessian: public SimulationStepBase_EnergyForceTorque{

    protected:

        uninitialized_cached_vector<tensor3> hessianTmp;

        void copyToTmp(cudaStream_t st){

                SimulationStepBase_EnergyForceTorque::copyToTmp(st);

                System::log<System::DEBUG>("[SimulationStep_EnergyForceTorqueHessian] (%s) Copying to tmp", name.c_str());

                int N = this->pd->getNumParticles();

                hessianTmp.resize(N);

                {
                        auto hessian = this->pd->getHessian(access::location::gpu, access::mode::read);
                        thrust::copy(thrust::cuda::par.on(st),
                                     hessian.begin(),
                                     hessian.end(),
                                     hessianTmp.begin());
                }

                cudaDeviceSynchronize();

        }

        void copyFromTmp(cudaStream_t st){

                SimulationStepBase_EnergyForceTorque::copyFromTmp(st);

                System::log<System::DEBUG>("[SimulationStep_EnergyForceTorqueHessian] (%s) Copying from tmp", name.c_str());

                int N = this->pd->getNumParticles();

                {
                        auto hessian = this->pd->getHessian(access::location::gpu, access::mode::write);
                        thrust::copy(thrust::cuda::par.on(st),
                                     hessianTmp.begin(),
                                     hessianTmp.end(),
                                     hessian.begin());
                }

                cudaDeviceSynchronize();

        }

        void setZero(cudaStream_t st){

                SimulationStepBase_EnergyForceTorque::setZero(st);

                System::log<System::DEBUG>("[SimulationStep_EnergyForceTorqueHessian] (%s) Setting to zero", name.c_str());

                {
                        auto hessian = this->pd->getHessian(access::location::gpu, access::mode::write);
                        thrust::fill(thrust::cuda::par.on(st),
                                     hessian.begin(),
                                     hessian.end(),
                                     tensor3(0.0));
                }

                cudaDeviceSynchronize();

        }

    public:

        SimulationStepBase_EnergyForceTorqueHessian(std::shared_ptr<ParticleGroup>              pg,
                                                    std::shared_ptr<IntegratorManager>  integrator,
                                                    std::shared_ptr<ForceField>                 ff,
                                                    DataEntry& data,
                                                    std::string name):SimulationStepBase_EnergyForceTorque(pg,integrator,ff,data,name){}

        ~SimulationStepBase_EnergyForceTorqueHessian(){}

        void tryApplyStep(cudaStream_t st,bool force=false){
                if(this->isStepApplied(force)){
                        this->copyToTmp(st);
                        SimulationStepBase_EnergyForceTorque::tryApplyStep(st,true);
                        this->copyFromTmp(st);
                }
        }

        void tryApplyStep(bool force=false){
                if(this->isStepApplied(force)){
                        this->copyToTmp(0);
                        SimulationStepBase_EnergyForceTorque::tryApplyStep(true);
                        this->copyFromTmp(0);
                }
        }
};

}}}

#endif
