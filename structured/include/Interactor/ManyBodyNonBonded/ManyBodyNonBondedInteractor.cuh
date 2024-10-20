#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "Interactor/ManyBodyNonBonded/ManyBodyNonBondedLoader.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Types.cuh"

namespace uammd {
namespace structured {
namespace Interactor {

    template<int THREADS_PER_BLOCK = 256>
    class ManyBodyNonBonded_ : public Interactor {
        protected:
            std::shared_ptr<GlobalData> gd;
            std::shared_ptr<ExtendedParticleData> pd;

            std::shared_ptr<InputEntryManager> manyBodyForceFieldInfo;

        public:
            ManyBodyNonBonded_(std::shared_ptr<GlobalData> gd,
                               std::shared_ptr<ParticleGroup> pg,
                               std::vector<std::string> path,
                               std::string name);

            virtual void sum(uammd::Interactor::Computables comp, cudaStream_t st) override;
    };

}}}
