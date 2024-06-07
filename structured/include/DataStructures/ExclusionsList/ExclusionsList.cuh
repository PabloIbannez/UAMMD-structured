// Header guards to prevent double inclusion of this header file
#ifndef __EXCLUSION_LIST__
#define __EXCLUSION_LIST__

namespace uammd{
namespace structured{

    // A struct representing a list of particles that are to be excluded for some operation
    struct particleExclusionList{

        // Two-dimensional pointer for exclusionStart, points to start of exclusion list for each particle
        int** exclusionStart;

        // Pointer to array that stores number of exclusions per particle
        int*  nExclusionPerParticle;

        // Indices for current particle and current exclusion
        int current_i;
        int current_nExc;

        // Pointers to start of exclusion list and buffer for current particle
        int* start;
        int* startBuffer;

        // Constructor for the particleExclusionList struct
        particleExclusionList(int** exclusionStart,
                              int*  nExclusionPerParticle):exclusionStart(exclusionStart),
                                                           nExclusionPerParticle(nExclusionPerParticle){}

        // Function to set the current particle
        inline __device__ void set(int i){
            current_i=i;
            current_nExc=nExclusionPerParticle[current_i];
            start = exclusionStart[current_i];

            startBuffer=start;
        }

        // Function to set the current particle and startBuffer, which is a pointer to shared memory
        // Shared memory in GPUs is a type of local memory that can be used for faster data access.
        // This is because shared memory has lower latency and higher bandwidth compared to global memory in GPUs.
        // Storing the exclusion list in shared memory improves performance as it speeds up subsequent accesses
        // to the exclusion list in the isPartExcluded function.
        // Also, if the exclusion lists for different particles have overlapping entries, storing them in shared memory
        // allows these common entries to be loaded once and then reused by multiple threads, saving on memory bandwidth.
        inline __device__ void set(int i,int* strB){
            current_i=i;
            current_nExc=nExclusionPerParticle[current_i];
            start = exclusionStart[current_i];
            startBuffer = strB;

            // Copying the exclusion list to the buffer in shared memory
            forj(0,current_nExc){
                startBuffer[j]=start[j];
            }
        }

        // Function to check if a particle is in the exclusion list
        inline __device__ bool isPartExcluded(int k){

            // Search the exclusion list for particle k
            forj(0,current_nExc){
                // If particle k is found, return true
                if(startBuffer[j]==k){return true;}
            }

            // If particle k is not in the exclusion list, return false
            return false;
        }

    };

    // Class representing the full list of particle exclusions
    class Exclusions{

        private:

            // Shared pointer to GlobalData
            std::shared_ptr<GlobalData>            gd;

            // Shared pointer to ExtendedParticleData
            std::shared_ptr<ExtendedParticleData>  pd;

            // Shared pointer to groupsList
            std::shared_ptr<groupsList> groupsList;

        public:

            // Constructor for Exclusions class
            Exclusions(std::shared_ptr<GlobalData>            gd,
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

            // Function to get the maximum number of exclusions
            int getMaxExclusions(){return groupsList->getMaxGroupSize();}

            // Function to get the shared memory size required for the exclusion list
            size_t getSharedSize(){return this->getMaxExclusions()*sizeof(int);}

            // Function to generate a particleExclusionList
            particleExclusionList getParticleExclusionList(){

                // Get the groupsListInfo from the groupsList
                auto groupsListInfo = groupsList->getGroupsListInfo(access::location::gpu);

                // Create and return a particleExclusionList with this information
                return particleExclusionList(groupsListInfo.listStart,
                                             groupsListInfo.n);
            }
    };

}}

// End of header guards
#endif
