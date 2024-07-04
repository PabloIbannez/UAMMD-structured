#pragma once

//TODO Check groups

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "Interactor/Patches/GenericPatchesPotentialLoader.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Types.cuh"

#include "utils/quaternion.cuh"

namespace uammd {
namespace structured {
namespace Interactor {

    template<int THREADS_PER_BLOCK = 256>
    class PatchyParticles_ : public Interactor {
        protected:
            std::shared_ptr<GlobalData> gd;
            std::shared_ptr<ExtendedParticleData> pd;

            std::shared_ptr<GlobalData> patchesGd;
            std::shared_ptr<ExtendedParticleData> patchesPd;

            std::shared_ptr<InputEntryManager> patchForceFieldInfo;

            std::map<std::string, std::shared_ptr<ParticleGroup>> patchesGroups;
            std::map<std::string, std::shared_ptr<VerletConditionalListSetBase>> patchesVConListSet;
            std::map<std::string, std::shared_ptr<Interactor>> patchInteractors;

            std::shared_ptr<groupsList> particleId2patchesId;

            void accumulate(uammd::Interactor::Computables comp, cudaStream_t st);

        public:
            PatchyParticles_(std::shared_ptr<GlobalData> gd,
                             std::shared_ptr<ParticleGroup> pg,
                             std::vector<std::string> path,
                             std::string name);

            void updatePatchyParticles(cudaStream_t st);

            std::shared_ptr<GlobalData> getPatchesGlobalData() { return patchesGd; }
            std::shared_ptr<ExtendedParticleData> getPatchesParticleData() { return patchesPd; }
            std::map<std::string, std::shared_ptr<Interactor>> getPatchesInteractors() { return patchInteractors; }

            InputEntryManager::entryInfo getPatchesInteractorInfo(std::string name);

            virtual void sum(uammd::Interactor::Computables comp, cudaStream_t st) override;
    };

    template<int THREADS_PER_BLOCK = 256>
    class DynamicallyBondedPatchyParticles_ : public PatchyParticles_<THREADS_PER_BLOCK> {
        private:
            using Base = PatchyParticles_<THREADS_PER_BLOCK>;
            real energyThreshold;

        public:
            DynamicallyBondedPatchyParticles_(std::shared_ptr<GlobalData> gd,
                    std::shared_ptr<ParticleGroup> pg,
                    std::vector<std::string> path,
                    std::string name);

            void updatePatchyParticles(cudaStream_t st);
            void sum(uammd::Interactor::Computables comp, cudaStream_t st) override;
    };

}}}
