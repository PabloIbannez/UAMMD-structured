#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Single/SingleInteractor.cuh"
#include "Interactor/Single/External/External.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "ExternalTabulated.cuh"

REGISTER_SINGLE_INTERACTOR(
    External,ExternalTabulated,
    uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::External::ExternalTabulated>
)
