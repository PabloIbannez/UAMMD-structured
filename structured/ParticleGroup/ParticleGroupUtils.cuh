#ifndef __GROUP_UTILS__
#define __GROUP_UTILS__

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


    std::map<int,std::shared_ptr<uammd::ParticleGroup>> BatchGroupLoader(std::shared_ptr<ParticleGroup> pg){

        std::shared_ptr<ParticleData> pd = pg->getParticleData();

        std::map<int,std::shared_ptr<uammd::ParticleGroup>> batGroup;

        std::set<int> batchList;
        {
            auto groupIndex = pg->getIndexIterator(access::location::cpu);

            auto batchId = pd->getBatchId(uammd::access::location::cpu,uammd::access::mode::read);

            fori(0,pg->getNumberParticles()){
                batchList.emplace(batchId[groupIndex[i]]);
            }
        }

        std::map<int,std::vector<int>> batchIdIds;
        {
            auto groupIndex = pg->getIndexIterator(access::location::cpu);

            auto batchId = pd->getBatchId(uammd::access::location::cpu,uammd::access::mode::read);
            auto id    = pd->getId(uammd::access::location::cpu,uammd::access::mode::read);

            for(int bId : batchList){
                fori(0,pg->getNumberParticles()){
                    if (bId == batchId[groupIndex[i]]){
                        batchIdIds[bId].push_back(id[groupIndex[i]]);
                    }
                }

            }
        }

        for(const int& b : batchList){

            std::string groupName = "batchId_"+std::to_string(b);

            auto pg = std::make_shared<uammd::ParticleGroup>(batchIdIds[b].begin(),batchIdIds[b].end(),
                                                             pd,
                                                             groupName);
            batGroup[b]=pg;

        }

        return batGroup;
    }

    std::vector<int> getBatchesInGroup(std::shared_ptr<ParticleGroup> pg){

        std::shared_ptr<ParticleData> pd = pg->getParticleData();

        std::set<int> batchList;
        {
            auto groupIndex = pg->getIndexIterator(access::location::cpu);
            auto batchId = pd->getBatchId(uammd::access::location::cpu,uammd::access::mode::read);

            fori(0,pg->getNumberParticles()){
                batchList.emplace(batchId[groupIndex[i]]);
            }
        }

        return std::vector<int>(batchList.begin(),batchList.end());
    }

    std::vector<int> getBatchesInParticleData(std::shared_ptr<ParticleData> pd){

        std::set<int> batchList;
        {
            auto batchId = pd->getBatchId(uammd::access::location::cpu,uammd::access::mode::read);

            fori(0,pd->getNumParticles()){
                batchList.emplace(batchId[i]);
            }
        }

        return std::vector<int>(batchList.begin(),batchList.end());
    }

    int BatchGroupNumber(std::shared_ptr<ParticleGroup> pg){
        return getBatchesInGroup(pg).size();
    }

    //Stop simulation if BatchGroupNumber is larger than n
    void BatchGroupNumberCheck(std::shared_ptr<ParticleGroup> pg, int n){
        int batchGroupNumber = BatchGroupNumber(pg);
        if( batchGroupNumber > n){
            System::log<System::CRITICAL>("[BatchGroupNumberCheck] The number of batches in the group"
                                          " is larger than the limit (%d > %d)." , batchGroupNumber, n);
        }
    }

    int BatchNumber(std::shared_ptr<ParticleData> pd){
        return getBatchesInParticleData(pd).size();
    }

    std::map<std::string,std::shared_ptr<ParticleGroup>> loadGroupsList(std::shared_ptr<ExtendedSystem>      sys,
                                                                        std::shared_ptr<GlobalData>           gd,
                                                                        std::shared_ptr<ExtendedParticleData> pd,
                                                                        std::vector<std::string> path){

        auto data = sys->getInput()->getDataEntry(path);

        std::string name = path.back();

        if(!((data.getType() == "Groups") and (data.getSubType() == "GroupsList"))){
            System::log<System::CRITICAL>("[GroupLoader] Expected entry \"%s\" to be a \"Groups\" entry type and \"GroupsList\" sub type."
                                          " But entry type is: %s, and sub type: %s",
                                          data.name.c_str(),
                                          data.getType().c_str(),
                                          data.getSubType().c_str());
        }

        std::vector<std::string> availSelectors = {"Ids","notIds",
                                                   "Types","notTypes",
                                                   "ModelIds","notModelIds",
                                                   "BatchIds","notBatchIds",
                                                   "ModelIdsBatchIds","notModelIdsBatchIds"};

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;

        auto groupsData = data.getDataMap();

        for(auto& gData : groupsData){

            std::string gType = gData.at("type");
            std::string gName = gData.at("name");

            System::log<System::MESSAGE>("[GroupLoader] (%s) Loading group \"%s\" of type \"%s\".",
                                         name.c_str(),
                                         gName.c_str(),
                                         gType.c_str());


            if        (gType == "Ids" or gType == "notIds"){
                std::vector<int> idsList = gData.at("selection");

                if(gType == "Ids"){
                    auto sel = selectors::ids(idsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                } else {
                    auto sel = selectors::notIds(idsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                }

            } else if (gType == "Types" or gType == "notTypes"){

                std::vector<std::string> selectedTypeNamesList = gData.at("selection");

                //Name to type id
                auto typesHandler = gd->getTypes();

                std::vector<int> typesList;
                for(std::string typeName : selectedTypeNamesList){
                    typesList.push_back(typesHandler->getTypeId(typeName));
                }

                if (gType == "Types"){
                    auto sel = selectors::types(typesList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                } else {
                    auto sel = selectors::notTypes(typesList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                }

            } else if (gType == "ModelIds" or gType == "notModelIds"){

                std::vector<int> modelIdsList = gData.at("selection");

                if(gType == "ModelIds"){
                    auto sel = selectors::modelIds(modelIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                } else {
                    auto sel = selectors::notModelIds(modelIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                }

            } else if (gType == "BatchIds" or gType == "notBatchIds"){

                std::vector<int> batchIdsList = gData.at("selection");

                if(gType == "BatchIds"){
                    auto sel = selectors::batchIds(batchIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                } else {
                    auto sel = selectors::notBatchIds(batchIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                }

            } else if (gType == "ModelIdsBatchIds" or gType == "notModelIdsBatchIds"){

                std::vector<std::vector<int>> mdlIdsBatchIdsList = gData.at("selection");

                std::vector<int> mdlIdsList;
                std::vector<int> batchIdsList;
                for(auto& mdlBatchIds : mdlIdsBatchIdsList){
                    mdlIdsList.push_back(mdlBatchIds[0]);
                    batchIdsList.push_back(mdlBatchIds[1]);
                }

                if(gType == "ModelIdsBatchIds"){
                    auto sel = selectors::modelIdsBatchIds(mdlIdsList,batchIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                } else {
                    auto sel = selectors::notModelIdsBatchIds(mdlIdsList,batchIdsList);
                    groups[gName] = std::make_shared<ParticleGroup>(sel,pd,gName);
                }

            } else {

                bool first = true;
                std::string message = "";
                for(std::string selType: availSelectors){
                    if (first){
                        first = false;
                    } else {
                        message+=", ";
                    }
                    message+=selType;
                }

                System::log<System::CRITICAL>("[GroupLoader] (%s) Selection \"%s\" is not available."
                                              " Available selections are: %s",name.c_str(),gType.c_str(),message.c_str());
            }



        }

        return groups;
    }

    std::map<std::string,std::shared_ptr<ParticleGroup>> loadGroupsListFromInputEntries(std::shared_ptr<ExtendedSystem>        sys,
                                                                                        std::shared_ptr<GlobalData>             gd,
                                                                                        std::shared_ptr<ExtendedParticleData>   pd,
                                                                                        std::shared_ptr<InputEntryManager> entries){

        std::map<std::string,std::shared_ptr<ParticleGroup>> groups;

        groups["All"] = std::make_shared<ParticleGroup>(pd,"All");

        for(auto& entry : entries->getEntriesInfo()){
            if(entry.second.entryType == "Groups" and entry.second.entrySubType == "GroupsList"){
                std::map<std::string,std::shared_ptr<ParticleGroup>> groupsTmp = loadGroupsList(sys,gd,pd,entry.second.path);

                //Check if some group is already present
                for(auto g : groupsTmp){
                    if(groups.count(g.first) == 0){
                        groups[g.first] = g.second;
                        entry.second.used = true;
                    } else {
                        System::log<System::CRITICAL>("[GroupsListFromInputEntries] Error loading groups,"
                                                      "group \"%s\" has already been added.",g.first.c_str());
                    }
                }
            }
        }

        //Print information about groups
        for(auto g : groups){
            System::log<System::MESSAGE>("[GroupsListFromInputEntries] Added group \"%s\" with %i particles.",
                                         g.first.c_str(),g.second->getNumberParticles());
        }

        return groups;
    }

    std::shared_ptr<ParticleGroup> getParticleGroupFromGroupsList(std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                                                                  DataEntry& data,
                                                                  std::string defaultGroup){
        //Load group
        std::string groupName;
        if(not data.isParameterAdded("group")){
            groupName = defaultGroup;
            System::log<System::MESSAGE>("[ParticleGroupFromGroupsList] No group name specified, using \"%s\".",groupName.c_str());
        } else {
            groupName = data.getParameter<std::string>("group");
        }
        System::log<System::MESSAGE>("[ParticleGroupFromGroupsList] Using group \"%s\".",groupName.c_str());

        //Check if the group exists
        if(groups.find(groupName) == groups.end()){
            System::log<System::CRITICAL>("[ParticleGroupFromGroupsList] Group \"%s\" not found.",groupName.c_str());
        }

        return groups[groupName];
    }
}}}

#endif
