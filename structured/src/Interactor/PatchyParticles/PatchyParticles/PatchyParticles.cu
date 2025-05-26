#pragma once
#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/PatchyParticles/PatchyParticlesInteractor.cuh"
#include "Interactor/InteractorFactory.cuh"

namespace uammd{
namespace structured{
namespace Interactor{
namespace PatchyParticles{

    using PatchyParticles = typename Interactor::PatchyParticles_<>;

}}}}

REGISTER_PATCHY_PARTICLE_INTERACTOR(
    PatchyParticles,PatchyParticles,
    uammd::structured::Interactor::PatchyParticles::PatchyParticles
)
