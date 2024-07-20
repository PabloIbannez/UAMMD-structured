#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/UnitsFactory.cuh"

namespace uammd{
namespace structured{
namespace UnitsLoader{

    bool isUnitsAvailable(std::shared_ptr<ExtendedSystem> sys,
                          std::vector<std::string>       path);

    std::shared_ptr<typename Units::UnitsHandler>
    loadUnits(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path);

}}}
