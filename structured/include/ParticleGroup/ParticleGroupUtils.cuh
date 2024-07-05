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

    namespace selectors{

        namespace selectors_ns{

            template<class sel>
            class selector_negation: public sel{
                public:
                    //inheriting constructor
                    using sel::sel;

                    bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                        return !sel::isSelected(particleIndex,pd);
                    }
            };

        }

        class ids{

            private:

                std::vector<int> idsList;

            public:

                ids(const std::vector<int>& idsList):idsList(idsList){};

                bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                    int id = (pd->getId(access::cpu, access::read).raw())[particleIndex];

                    bool found = (std::find(idsList.begin(),idsList.end(),id) != idsList.end());

                    return found;
                }
        };

        using notIds = selectors_ns::selector_negation<ids>;

        class types{

            private:

                std::vector<int> typesList;

            public:

                types(const std::vector<int>& typesList):typesList(typesList){};

                bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                    int type = int((pd->getPos(access::cpu, access::read).raw())[particleIndex].w);

                    bool found = (std::find(typesList.begin(),typesList.end(),type) != typesList.end());

                    return found;
                }
        };

        using notTypes = selectors_ns::selector_negation<types>;

        class modelIds{

            private:

                std::vector<int> modelIdsList;

            public:

                modelIds(const std::vector<int>& modelIdsList):modelIdsList(modelIdsList){};

                bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                    int mdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];

                    bool found = (std::find(modelIdsList.begin(),modelIdsList.end(),mdlId) != modelIdsList.end());

                    return found;
                }
        };

        using notModelIds = selectors_ns::selector_negation<modelIds>;

        class batchIds{

            private:

                std::vector<int> batchIdsList;

            public:

                batchIds(const std::vector<int>& batchIdsList):batchIdsList(batchIdsList){};

                bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                    int batchId = (pd->getBatchId(access::cpu, access::read).raw())[particleIndex];

                    bool found = (std::find(batchIdsList.begin(),batchIdsList.end(),batchId) != batchIdsList.end());

                    return found;
                }
        };

        using notBatchIds = selectors_ns::selector_negation<batchIds>;

        class modelIdsBatchIds{

            private:

                int n;

                std::vector<int> modelIdsList;
                std::vector<int> batchIdsList;

            public:

                modelIdsBatchIds(const std::vector<int>& modelIdsList,
                               const std::vector<int>& batchIdsList):modelIdsList(modelIdsList),
                                                                   batchIdsList(batchIdsList){
                    //Check if the sizes of the two vectors are the same
                    if(modelIdsList.size() != batchIdsList.size()){
                        System::log<System::CRITICAL>("[Selector ModelIdsBatchIds] The size of the modelIdsList and batchIdsList must be the same.");
                    }

                    n = modelIdsList.size();
                };

                bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                    int mdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];
                    int batchId = (pd->getBatchId(access::cpu, access::read).raw())[particleIndex];

                    bool found = false;
                    for(int i=0;i<n and !found;i++){
                        if(mdlId == modelIdsList[i] and batchId == batchIdsList[i]){
                            found = true;
                        }
                    }

                    return found;
                }
        };

        using notModelIdsBatchIds = selectors_ns::selector_negation<modelIdsBatchIds>;

    }

    template<class selector>
    std::shared_ptr<uammd::ParticleGroup> createSubGroup(std::shared_ptr<uammd::ParticleGroup> pg,
                                                         selector sel,
                                                         std::string newGroupName){

        std::vector<int> selectedParticles;

        auto pd = pg->getParticleData();

        {
            auto groupIndex = pg->getIndexIterator(access::location::cpu);
            auto ids = pd->getId(uammd::access::location::cpu,uammd::access::mode::read);

            fori(0,pg->getNumberParticles()){
                if(sel.isSelected(groupIndex[i],pd)){
                    selectedParticles.push_back(ids[groupIndex[i]]);
                }
            }
        }

        auto subGroup = std::make_shared<uammd::ParticleGroup>(selectedParticles.begin(),selectedParticles.end(),
                                                               pd,
                                                               newGroupName);

        return subGroup;
    }

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
