#include "DataStructures/ExclusionsList/ExclusionsList.cuh"

namespace uammd{
namespace structured{

    Exclusions::Exclusions(std::shared_ptr<GlobalData>            gd,
                           std::shared_ptr<ExtendedParticleData>  pd,
                           DataEntry& dataEntry):gd(gd),pd(pd){

        System::log<System::MESSAGE>("[Exclusions] Initialized");

        // Get total number of particles
        int N = pd->getNumParticles();

        // Load the groupsList from the dataEntry
        groupsList = dataEntry.getDataAsGroupsList("id","id_list",0,N-1,true);
        // Note: for a system with N particles, the id go from 0 to N-1
        // Note: true flag is for checking that relations are symmetric

        // Check if the number of groups is the same as the number of particles
        if(N!=groupsList->getNGroups()){
            System::log<System::CRITICAL>("[Exclusions] The number of entries"
                                          " in the exclusions list differs "
                                          "from the number of particles");
        }

    }

    // Function to generate a particleExclusionList
    particleExclusionList Exclusions::getParticleExclusionList(){

        // Get the groupsListInfo from the groupsList
        auto groupsListInfo = groupsList->getGroupsListInfo(access::location::gpu);

        // Create and return a particleExclusionList with this information
        return particleExclusionList(groupsListInfo.listStart,
                                     groupsListInfo.n);
    }

}}
