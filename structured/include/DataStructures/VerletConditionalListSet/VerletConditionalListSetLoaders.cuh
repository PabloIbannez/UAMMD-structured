#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetFactory.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

namespace uammd{
namespace structured{
namespace VerletConditionalListSetLoaders{

    std::shared_ptr<uammd::structured::VerletConditionalListSetBase>
    loadVerletConditionalListSet(std::shared_ptr<ExtendedSystem> sys,
                                 std::shared_ptr<GlobalData>     gd,
                                 std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                 std::vector<std::string>       path);
}}}
