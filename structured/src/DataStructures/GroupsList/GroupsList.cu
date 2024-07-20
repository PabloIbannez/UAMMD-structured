#include "DataStructures/GroupsList/GroupsList.cuh"

namespace uammd{
namespace structured{

    void groupsList::init(const std::vector<std::vector<int>>& groups, bool sym){

         // Generate a sequence of integer IDs for the groups
         std::vector<int> id(groups.size());
         std::iota(id.begin(),id.end(),0);

         // Call the overloaded init function with the group IDs and groups
         init(id,groups,sym);
    }

    // Another 'init' function that takes a list of IDs and a nested vector of groups.
    void groupsList::init(std::vector<int> id,
                          std::vector<std::vector<int>> groups, bool sym){

         // Check ids are unique
         if(std::set<int>(id.begin(),id.end()).size() != id.size()){
             System::log<System::CRITICAL>("[GroupsList] The IDs are not unique");
         }

         // If the number of IDs does not match the number of groups, log a critical error.
         if(id.size() != groups.size()){
             System::log<System::CRITICAL>("[GroupsList] The size of id (%i) "
                     "does not match the size of groups (%i)",
                     id.size(),groups.size());
         }

         // id is sort in ascending order, and the groups are sorted for keeping the previous id-group relation.
         std::vector<std::pair<int,std::vector<int>>> idGroups;
         for(uint i = 0; i < id.size(); ++i){
             idGroups.push_back(std::make_pair(id[i],groups[i]));
         }
         std::sort(idGroups.begin(),idGroups.end());
         for(uint i = 0; i < id.size(); ++i){
             id[i]     = idGroups[i].first;
             groups[i] = idGroups[i].second;
         }

         // Ensures the IDs are a sequence from 0 to n-1, where n is the number of groups.
         // If not, logs a message and throws a domain error.
         {

             auto minmax = std::minmax_element(id.begin(), id.end());

             auto min = minmax.first;
             auto max = minmax.second;

             std::vector<int> toCheck(*max+1);
             std::iota(toCheck.begin(),toCheck.end(),*min);

             try{
                 for(uint i=0;i<id.size();i++){
                     if(toCheck[i]!=id[i]){
                         System::log<System::ERROR>("[GroupsList] The id %i"
                                 " is not in the current id vector",
                                 toCheck[i]);
                         throw std::domain_error(std::to_string(toCheck[i]));
                     }
                 }
             } catch (...){
                 System::log<System::CRITICAL>("[GroupsList] Id vector error");
             }
         }

         // Counts the total number of groups, total number of relations, and the maximum group size.

         nGroups      = id.size();
         nRelations   = 0;
         maxGroupSize = 0;

         for(const std::vector<int>& grp : groups){
             nRelations+=grp.size();
             maxGroupSize = std::max(maxGroupSize,uint(grp.size()));
         }

         System::log<System::MESSAGE>("[GroupsList] Number of groups: %i, number of relations: %i, "
                 "max group size: %i",
                 nGroups,nRelations,maxGroupSize);

         if(sym){
             // We check if relations are symmetric, i.e. if a is in the group of b, then b is in the group of a.
             // If not, logs a message and throws a domain error.

             for(uint i = 0; i < id.size(); ++i){
                 for(uint j = 0; j < groups[i].size(); ++j){
                     int a = id[i];
                     int b = groups[i][j];
                     if(std::find(groups[b].begin(),groups[b].end(),a) == groups[b].end()){
                         System::log<System::MESSAGE>("[GroupsList] The relation (%i,%i) is not symmetric",
                                 a,b);
                         throw std::domain_error(std::to_string(a)+","+std::to_string(b));
                     }
                 }
             }
         }

         // Constructs the host vectors for list, listStart and n, and fills in the data.

         h_list.resize(nRelations);
         h_listStart.resize(nGroups,nullptr);
         h_n.resize(nGroups,0);

         int relationCounter = 0;
         for(uint i=0;i<nGroups;i++){
             std::vector<int> currentGroup = groups[i];
             for(uint j=0;j<currentGroup.size();j++){
                 h_list[relationCounter] = currentGroup[j];
                 relationCounter++;
             }
             h_n[i] = currentGroup.size();
         }

         int offSet=0;
         for(uint i=0;i<nGroups;i++){
             int currentId = id[i];
             if(h_n[currentId]!=0){
                 h_listStart[currentId]=h_list.data()+offSet;
                 offSet=offSet+h_n[currentId];
             }
         }

         //for(uint i=0;i<nGroups;i++){
         //    int currentId = id[i];
         //    std::cout << currentId << ": ";
         //    for(uint n=0;n<h_n[currentId];n++){
         //        int value = h_listStart[currentId][n];
         //        std::cout << value << " ";
         //    }
         //    std::cout << std::endl;
         //}
         //std::cin.get();

         // Copies the host vectors to the device vectors.

         list = h_list;
         n    = h_n;

         thrust::host_vector<int*> tmp = h_listStart;
         offSet=0;
         for(uint i=0;i<nGroups;i++){
             if(h_n[i]!=0){
                 h_listStart[i]=thrust::raw_pointer_cast(list.data())+offSet;
                 offSet=offSet+h_n[i];
             }
         }
         listStart = h_listStart;
         h_listStart = tmp;
     }

     // Constructor that initializes the group list from a nested vector of groups.
    groupsList::groupsList(std::vector<std::vector<int>> groups,bool sym){ init(groups,sym); }

     // Constructor that initializes the group list from a nested vector of groups.
    groupsList::groupsList(std::vector<int> id,
                           std::vector<std::vector<int>> groups, bool sym){ init(id,groups,sym); }

    groupsList::groupsListInfo groupsList::getGroupsListInfo(access::location loc){

        groupsList::groupsListInfo gListInfo;

        if        (loc==access::location::cpu){
            gListInfo.list      = thrust::raw_pointer_cast(h_list.data());
            gListInfo.listStart = thrust::raw_pointer_cast(h_listStart.data());
            gListInfo.n         = thrust::raw_pointer_cast(h_n.data());
        } else if (loc==access::location::gpu){
            gListInfo.list      = thrust::raw_pointer_cast(list.data());
            gListInfo.listStart = thrust::raw_pointer_cast(listStart.data());
            gListInfo.n         = thrust::raw_pointer_cast(n.data());
        } else {
            System::log<System::CRITICAL>("[GroupsList] Selected access::location not avaible");
        }

        gListInfo.nGroups      = nGroups;
        gListInfo.nRelations   = nRelations;
        gListInfo.maxGroupSize = maxGroupSize;

        return gListInfo;
   }
}}
