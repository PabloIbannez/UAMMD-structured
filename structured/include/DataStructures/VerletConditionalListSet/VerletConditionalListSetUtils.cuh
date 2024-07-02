#pragma once

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetLoaders.cuh"

namespace uammd{
namespace structured{
namespace VerletConditionalListSetUtils{

    std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>
    loadVerletConditionalListSetsFromInputEntries(std::shared_ptr<ExtendedSystem>        sys,
                                                  std::shared_ptr<GlobalData>             gd,
                                                  std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                                  std::shared_ptr<InputEntryManager> entries);

    std::shared_ptr<VerletConditionalListSetBase>
    getNeighbourListFromNeighbourListsList(std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                                           DataEntry& data);

}}}
