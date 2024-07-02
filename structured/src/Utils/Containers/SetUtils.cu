#include "Utils/Containers/SetUtils.cuh"

namespace uammd{
namespace structured{

std::vector<int> setsIntersection(const std::vector<std::vector<std::vector<int>>>& idSets){

    std::vector<int> allIds;

    std::vector<int> intersection;

    for(const std::vector<std::vector<int>>& idSet : idSets){
        for(const std::vector<int>& idSubSet : idSet){
            std::vector<int> setBuffer = idSubSet;

            std::sort(setBuffer.begin(),setBuffer.end());
            std::sort(allIds.begin(),allIds.end());

            std::set_intersection(setBuffer.begin(),setBuffer.end(),
                    allIds.begin(),allIds.end(),
                    std::back_inserter(intersection));

            allIds.insert(allIds.end(),idSubSet.begin(),idSubSet.end());
        }
    }

    return intersection;
}

}}
