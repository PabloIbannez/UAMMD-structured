#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/Ensemble/EnsembleHandler.cuh"
#include "GlobalData/Ensemble/EnsembleFactory.cuh"

namespace uammd{
namespace structured{
namespace EnsembleLoader{

    bool isEnsembleAvailable(std::shared_ptr<ExtendedSystem> sys,
                             std::vector<std::string>       path);

    std::shared_ptr<typename Ensemble::EnsembleHandler>
    loadEnsemble(std::shared_ptr<ExtendedSystem> sys,
                 std::vector<std::string>       path);

}}}
