#include "SimulationStep/SimulationStep.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{

    bool SimulationStepBase::isStepApplied(bool force){
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
	return false;
    }

    SimulationStepBase::SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
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

    SimulationStepBase::SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
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
    SimulationStepBase::SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
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

    void SimulationStepBase::tryInit( cudaStream_t st){
        if(!initialized){
            this->init(st);
            initialized=true;
        }
    };

    void SimulationStepBase::tryInit(){
        cudaDeviceSynchronize();
        tryInit(0);
        cudaDeviceSynchronize();
    }


    void SimulationStepBase::tryApplyStep(cudaStream_t st,bool force){
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

    void SimulationStepBase::tryApplyStep(bool force){
        cudaDeviceSynchronize();
        tryApplyStep(0,force);
        cudaDeviceSynchronize();
    }

    //SimulationStepBase_EnergyForceTorque


    void SimulationStepBase_EnergyForceTorque::copyToTmp(cudaStream_t st){

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

    void SimulationStepBase_EnergyForceTorque::copyFromTmp(cudaStream_t st){

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

    void SimulationStepBase_EnergyForceTorque::setZero(cudaStream_t st){

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

    SimulationStepBase_EnergyForceTorque::SimulationStepBase_EnergyForceTorque
    (std::shared_ptr<ParticleGroup>              pg,
     std::shared_ptr<IntegratorManager>  integrator,
     std::shared_ptr<ForceField>                 ff,
     DataEntry& data,
     std::string name):SimulationStepBase(pg,integrator,ff,data,name){}

    void SimulationStepBase_EnergyForceTorque::tryApplyStep(cudaStream_t st,bool force){
        if(this->isStepApplied(force)){
            this->copyToTmp(st);
            SimulationStepBase::tryApplyStep(st,true);
            this->copyFromTmp(st);
        }
    }

    void SimulationStepBase_EnergyForceTorque::tryApplyStep(bool force){
        if(this->isStepApplied(force)){
            this->copyToTmp(0);
            SimulationStepBase::tryApplyStep(true);
            this->copyFromTmp(0);
        }
    }

    //SimulationStepBase_EnergyForceTorqueHessian


    void SimulationStepBase_EnergyForceTorqueHessian::copyToTmp(cudaStream_t st){

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

    void SimulationStepBase_EnergyForceTorqueHessian::copyFromTmp(cudaStream_t st){

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

    void SimulationStepBase_EnergyForceTorqueHessian::setZero(cudaStream_t st){

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

    SimulationStepBase_EnergyForceTorqueHessian::SimulationStepBase_EnergyForceTorqueHessian
    (std::shared_ptr<ParticleGroup>              pg,
     std::shared_ptr<IntegratorManager>  integrator,
     std::shared_ptr<ForceField>                 ff,
     DataEntry& data,
     std::string name):SimulationStepBase_EnergyForceTorque(pg,integrator,ff,data,name){}

    void SimulationStepBase_EnergyForceTorqueHessian::tryApplyStep(cudaStream_t st,bool force){
        if(this->isStepApplied(force)){
            this->copyToTmp(st);
            SimulationStepBase_EnergyForceTorque::tryApplyStep(st,true);
            this->copyFromTmp(st);
        }
    }

    void SimulationStepBase_EnergyForceTorqueHessian::tryApplyStep(bool force){
        if(this->isStepApplied(force)){
            this->copyToTmp(0);
            SimulationStepBase_EnergyForceTorque::tryApplyStep(true);
            this->copyFromTmp(0);
        }
    }

}}}
