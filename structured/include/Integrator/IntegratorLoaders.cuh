#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"

namespace uammd{
namespace structured{
namespace IntegratorLoader{

    bool isIntegratorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path);


    std::shared_ptr<typename uammd::Integrator>
    loadIntegrator(std::shared_ptr<ExtendedSystem> sys,
                   std::shared_ptr<GlobalData>     gd,
                   std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                   std::vector<std::string>       path);

}}}
