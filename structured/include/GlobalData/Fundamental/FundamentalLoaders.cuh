#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/FundamentalFactory.cuh"

namespace uammd{
namespace structured{
namespace FundamentalLoader{

    bool isFundamentalAvailable(std::shared_ptr<ExtendedSystem> sys,
                                std::vector<std::string>       path);

    std::shared_ptr<typename Fundamental::FundamentalHandler>
    loadFundamental(std::shared_ptr<ExtendedSystem> sys,
                    std::vector<std::string>       path);

}}}
