#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "Interactor/Patches/PatchesFactory.cuh"

#include <memory>
#include <vector>
#include <string>
#include <map>

namespace uammd{
namespace structured{
namespace Interactor{

    bool isPatchInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                                    std::vector<std::string>       path);

    std::shared_ptr<typename uammd::Interactor>
    loadPatchInteractor(std::shared_ptr<ExtendedSystem> sys,
                        std::shared_ptr<GlobalData>  gd,std::shared_ptr<ParticleGroup>  pg,
                        std::shared_ptr<GlobalData>  patchesGd,std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                        std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                        std::vector<std::string>       path);

}}}
