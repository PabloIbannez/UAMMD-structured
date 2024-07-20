#include "ParticleGroup/ParticleGroupUtils.cuh"

namespace uammd{
namespace structured{
namespace GroupUtils{

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
