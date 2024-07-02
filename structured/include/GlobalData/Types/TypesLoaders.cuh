#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/Types/TypesHandler.cuh"
#include "GlobalData/Types/TypesFactory.cuh"

namespace uammd{
namespace structured{
namespace TypesLoader{

    bool isTypesAvailable(std::shared_ptr<ExtendedSystem> sys,
                          std::vector<std::string>       path);

    std::shared_ptr<typename Types::TypesHandler>
    loadTypes(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path);

}}}
