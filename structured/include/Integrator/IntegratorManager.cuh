#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorUtils.cuh"
#include "Integrator/IntegratorLoaders.cuh"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace uammd{
namespace structured{

class IntegratorManager{

    private:

        std::shared_ptr<ExtendedSystem>       sys;
        std::shared_ptr<GlobalData>            gd;
        std::shared_ptr<ExtendedParticleData>  pd;

        std::vector<std::string> path;

        //////////////////////////////////////////

        std::shared_ptr<InputEntryManager>   integratorsInfo;

        //////////////////////////////////////////

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;
        std::map<std::string,std::shared_ptr<uammd::Integrator>> integrators;

        struct stepsInfo{
            std::string name;
            uint order;
            ullint steps;
        };

        std::map<std::string,stepsInfo> integratorSteps;

        //////////////////////////////////////////

        void loadGroups();
        void loadSchedule();

        void loadIntegrators();

    public:

        IntegratorManager(std::shared_ptr<ExtendedSystem>       sys,
                          std::shared_ptr<GlobalData>            gd,
                          std::shared_ptr<ExtendedParticleData>  pd,
                          std::vector<std::string> path);

        IntegratorManager(std::shared_ptr<ExtendedSystem>       sys,
                          std::shared_ptr<GlobalData>            gd,
                          std::shared_ptr<ExtendedParticleData>  pd);

        //Add a new interactor to the system
        void addIntegrator(std::shared_ptr<uammd::Integrator> integrator,
                           std::string name, int steps);

        std::map<std::string,std::shared_ptr<uammd::Integrator>>& getIntegrators();
        std::shared_ptr<uammd::Integrator> getIntegrator(std::string name);

        std::vector<stepsInfo> getSortedIntegratorSteps();
};

}}
