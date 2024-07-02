#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "GlobalData/GlobalData.cuh"

#include "Integrator/IntegratorBase.cuh"

namespace uammd{
namespace structured{
namespace IntegratorUtils{

    void generateVelocity(shared_ptr<ParticleGroup> pg,
                          real kBT,
                          uint seed,
                          cudaStream_t stream);

    void loadFrictionConstant(shared_ptr<ParticleGroup> pg,
                              real frictionConstant);

    void loadMobility(shared_ptr<ParticleGroup> pg,
                      real viscosity);

}}}
