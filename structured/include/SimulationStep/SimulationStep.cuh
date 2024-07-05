#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "ParticleGroup/ParticleGroupUtils.cuh"
#include "Integrator/IntegratorManager.cuh"
#include "ForceFields/ForceFields.cuh"

#include "Topology/Topology.cuh"

#include "utils/container.h"

#include "Utils/Backup/Backup.cuh"

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

        bool isStepApplied(bool force=false);
    public:

        SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           std::string name);

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
                           std::string name);

        //Construct from DataEntry
        SimulationStepBase(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           DataEntry& data,
                           std::string name);

        std::string getName(){return name;}

        ullint getStartStep()   {return startStep;}
        ullint getEndStep()     {return endStep;}
        ullint getIntervalStep(){return intervalStep;}

        ullint getLastStepApplied(){return lastStepApplied;}

        void tryInit( cudaStream_t st);
        void tryInit();

        void tryApplyStep(cudaStream_t st,bool force=false);
        void tryApplyStep(bool force=false);
};

class SimulationStepBase_EnergyForceTorque: public SimulationStepBase{

    protected:

        uninitialized_cached_vector<real>  energyTmp;
        uninitialized_cached_vector<real4> forceTmp;
        uninitialized_cached_vector<real4> torqueTmp;

        void copyToTmp(cudaStream_t st);
        void copyFromTmp(cudaStream_t st);
        void setZero(cudaStream_t st);

    public:

        SimulationStepBase_EnergyForceTorque(std::shared_ptr<ParticleGroup>              pg,
                std::shared_ptr<IntegratorManager>  integrator,
                std::shared_ptr<ForceField>                 ff,
                DataEntry& data,
                std::string name);

        void tryApplyStep(cudaStream_t st,bool force=false);

        void tryApplyStep(bool force=false);

};

class SimulationStepBase_EnergyForceTorqueHessian: public SimulationStepBase_EnergyForceTorque{

    protected:

        uninitialized_cached_vector<tensor3> hessianTmp;

        void copyToTmp(cudaStream_t st);
        void copyFromTmp(cudaStream_t st);
        void setZero(cudaStream_t st);

    public:

        SimulationStepBase_EnergyForceTorqueHessian(std::shared_ptr<ParticleGroup> pg,
                                                    std::shared_ptr<IntegratorManager>  integrator,
                                                    std::shared_ptr<ForceField>                 ff,
                                                    DataEntry& data,
                                                    std::string name);

        void tryApplyStep(cudaStream_t st,bool force=false);
        void tryApplyStep(bool force=false);

};

}}}
