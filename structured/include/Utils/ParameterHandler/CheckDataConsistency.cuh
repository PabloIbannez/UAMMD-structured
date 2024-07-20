#pragma once

#include "System/ExtendedSystem.cuh"
#include "ParticleData/ParticleGroup.cuh"

namespace uammd{
namespace structured{
namespace Batching{

    bool checkDataConsistency(DataEntry& data);

    bool checkParticleGroupDataConsistency(std::shared_ptr<ParticleGroup> pg,
                                           DataEntry& data);

}}}
