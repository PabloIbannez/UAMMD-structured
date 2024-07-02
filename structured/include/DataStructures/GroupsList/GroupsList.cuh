#pragma once

#include "uammd.cuh"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <numeric>
#include <vector>
#include <algorithm>

namespace uammd{
namespace structured{

    // The 'groupsList' class, to manage groups of related integer items.
    class groupsList{

           // Host vectors that store the actual list, the start indices of each group in the list, and the number of elements in each group.
           thrust::host_vector<int>  h_list;
           thrust::host_vector<int*> h_listStart;
           thrust::host_vector<int>  h_n;

           // Similar to host vectors, but these are for device.
           thrust::device_vector<int>  list;
           thrust::device_vector<int*> listStart;
           thrust::device_vector<int>  n;

           uint nGroups;      // Total number of groups.
           uint nRelations;   // Total number of elements across all groups.
           uint maxGroupSize; // The size of the largest group.

           // Private function for initializing the group list.
           void init(const std::vector<std::vector<int>>& groups, bool sym);

           // Another 'init' function that takes a list of IDs and a nested vector of groups.
           void init(std::vector<int> id,
                     std::vector<std::vector<int>> groups, bool sym);

        public:

            // A struct to hold all necessary group list information.
            struct groupsListInfo{

                int*   list;
                int**  listStart;
                int*   n;

                int nGroups;
                int nRelations;
                int maxGroupSize;
            };

            // Constructor that initializes the group list from a nested vector of groups.
            groupsList(std::vector<std::vector<int>> groups,bool sym = false);

            // Constructor that initializes the group list from a nested vector of groups.
            groupsList(std::vector<int> id,
                       std::vector<std::vector<int>> groups, bool sym = false);

            //Note: vectors are passed by value, so the original vectors are not modified.

            // Getter functions for the total number of groups, total number of relations, and the size of the largest group.
            int getNGroups()     {return nGroups;}
            int getNRelations()  {return nRelations;}
            int getMaxGroupSize(){return maxGroupSize;}

            // Function to get a 'groupsListInfo' object.
            // Depending on the location specified, it sets the list, listStart, and n pointers to point to either host or device memory.
            // Also sets the number of groups, number of relations, and maximum group size.
            groupsListInfo getGroupsListInfo(access::location loc);

    };

}}
