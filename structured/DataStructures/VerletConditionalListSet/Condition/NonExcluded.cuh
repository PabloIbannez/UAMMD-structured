#ifndef __CONDITION_NEXCL__
#define __CONDITION_NEXCL__

namespace uammd{
namespace structured{
namespace conditions{

    class nonExcluded : public excludedConditionBase{

        public:

            /////////////////////////////

            static const int condNum = 1;
            enum condIndex {NONEXCL=0};

            /////////////////////////////

            nonExcluded(std::shared_ptr<GlobalData>            gd,
                        std::shared_ptr<ExtendedParticleData>  pd,
                        DataEntry& dataEntry):excludedConditionBase(gd,pd,dataEntry){

                System::log<System::MESSAGE>("[Condition] Condition \"nonExcluded\" initialized");

            };

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="nonExcluded") {return NONEXCL;}

                System::log<System::CRITICAL>("[Condition] Requested a condition"
                                              " that is not present, %s",
                                              conditionName.c_str());

                return -1;

            }

            ////////////////////////////////////////////

            struct conditionChecker{

                int*  id;

                particleExclusionList prtExclList;
                int maxExclusions;

                inline __device__ void set(const int& i,
                                           const int& offsetBufferIndex,
                                           const char* sharedBuffer){
                    prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
                }

                inline __device__ void check(const int& i,const int& j,
                                             bool cond[condNum]){

                    for(int c=0;c<condNum;c++){
                        cond[c] = false;
                    }

                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[NONEXCL] = true;
                    }
                }

                conditionChecker(int* id,
                                 particleExclusionList prtExclList,
                                 int maxExclusions):id(id),
                                                    prtExclList(prtExclList),
                                                    maxExclusions(maxExclusions){}

            };

            conditionChecker getConditionChecker(){

                auto id      = pd->getId(access::location::gpu, access::mode::read);

                return conditionChecker(id.raw(),
                                        exclusionList->getParticleExclusionList(),
                                        exclusionList->getMaxExclusions());
            }

    };

}}}

#endif
