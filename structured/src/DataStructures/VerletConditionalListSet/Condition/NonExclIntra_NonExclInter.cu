#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSet.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetFactory.cuh"
#include "DataStructures/VerletConditionalListSet/Condition/Condition.cuh"

namespace uammd{
namespace structured{
namespace conditions{

    class nonExclIntra_nonExclInter : public excludedConditionBase {

        public:

            /////////////////////////////

            static const int condNum = 2;
            enum condIndex {INTRA=0,INTER=1};

            /////////////////////////////

            nonExclIntra_nonExclInter(std::shared_ptr<GlobalData>            gd,
                                      std::shared_ptr<ExtendedParticleData>  pd,
                                      DataEntry& dataEntry):excludedConditionBase(gd,pd,dataEntry){

                System::log<System::MESSAGE>("[Condition] Condition \"nonExclIntra_nonExclInter\" initialized");

            }

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="intra")  {return INTRA;}
                if     (conditionName=="inter")  {return INTER;}

                System::log<System::CRITICAL>("[Condition] Requested a condition"
                                              " that is not present, %s",
                                              conditionName.c_str());

                return -1;

            }

            ////////////////////////////////////////////

            struct conditionChecker{

                int*  id;
                int*  modelId;

                particleExclusionList prtExclList;
                int maxExclusions;

                inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                    prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
                }

                inline __device__ void check(const int& i,const int& j,bool cond[condNum]){

                    for(int c=0;c<condNum;c++){
                        cond[c] = false;
                    }

                    if(!prtExclList.isPartExcluded(id[j])){
                        if(modelId[i] == modelId[j]){
                            cond[INTRA]=true;
                        } else {
                            cond[INTER]=true;
                        }
                    }
                }

                conditionChecker(int* id,int* modelId,
                                 particleExclusionList prtExclList,
                                 int maxExclusions):
                                 id(id),modelId(modelId),
                                 prtExclList(prtExclList),
                                 maxExclusions(maxExclusions){}

            };


            conditionChecker getConditionChecker(){

                auto id      = pd->getId(access::location::gpu, access::mode::read);
                auto model   = pd->getModelId(access::location::gpu, access::mode::read);

                return conditionChecker(id.raw(),model.raw(),
                                        exclusionList->getParticleExclusionList(),
                                        exclusionList->getMaxExclusions());
            }
    };

}}}

REGISTER_VERLET_CONDITIONAL_LIST_SET(
    nonExclIntra_nonExclInter,
    uammd::structured::VerletConditionalListSet<uammd::structured::conditions::nonExclIntra_nonExclInter>
)
