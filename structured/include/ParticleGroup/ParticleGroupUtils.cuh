#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include <string>
#include <map>
#include <memory>

namespace uammd{
namespace structured{
namespace GroupUtils{

    template<class selector>
    std::shared_ptr<uammd::ParticleGroup> createSubGroup(std::shared_ptr<uammd::ParticleGroup> pg,
                                                         selector sel,
                                                         std::string newGroupName);

    std::map<int,std::shared_ptr<uammd::ParticleGroup>> BatchGroupLoader(std::shared_ptr<ParticleGroup> pg);

    std::vector<int> getBatchesInGroup(std::shared_ptr<ParticleGroup> pg);
    std::vector<int> getBatchesInParticleData(std::shared_ptr<ParticleData> pd);

    int BatchGroupNumber(std::shared_ptr<ParticleGroup> pg);

    //Stop simulation if BatchGroupNumber is larger than n
    void BatchGroupNumberCheck(std::shared_ptr<ParticleGroup> pg, int n);

    int BatchNumber(std::shared_ptr<ParticleData> pd);

    std::map<std::string,std::shared_ptr<ParticleGroup>> loadGroupsList(std::shared_ptr<ExtendedSystem>      sys,
                                                                        std::shared_ptr<GlobalData>           gd,
                                                                        std::shared_ptr<ExtendedParticleData> pd,
                                                                        std::vector<std::string> path);

    std::map<std::string,std::shared_ptr<ParticleGroup>> loadGroupsListFromInputEntries(std::shared_ptr<ExtendedSystem>        sys,
                                                                                        std::shared_ptr<GlobalData>             gd,
                                                                                        std::shared_ptr<ExtendedParticleData>   pd,
                                                                                        std::shared_ptr<InputEntryManager> entries);

    std::shared_ptr<ParticleGroup> getParticleGroupFromGroupsList(std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                                                                  DataEntry& data,
                                                                  std::string defaultGroup);

}}}
