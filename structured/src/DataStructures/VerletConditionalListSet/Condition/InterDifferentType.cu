#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSet.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetFactory.cuh"

namespace uammd{
namespace structured{
namespace conditions{

    class interDifferentType{

        private:

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;

        public:

            /////////////////////////////

            static const int condNum = 1;
            enum condIndex {INTER=0};

            /////////////////////////////


            interDifferentType(std::shared_ptr<GlobalData>            gd,
                               std::shared_ptr<ExtendedParticleData>  pd,
                               DataEntry& dataEntry):gd(gd),pd(pd){

                    System::log<System::MESSAGE>("[Condition] Condition \"interDifferentType\" initialized");
            }

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="inter") {return INTER;}

                System::log<System::CRITICAL>("[Condition] Requested a condition"
                                              " that is not present, %s",
                                              conditionName.c_str());
                return -1;
            }

            size_t getSharedSize(){
                return 0;
            }

            ////////////////////////////////////////////

            struct conditionChecker{

                real4*  pos;
                int*    modelId;

                inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){}

                inline __device__ void check(const int& i,const int& j,bool cond[condNum]){

                    for(int c=0;c<condNum;c++){
                        cond[c] = false;
                    }

                    if(modelId[i] != modelId[j] and int(pos[i].w) != int(pos[j].w)){
                        cond[INTER]=true;
                    }
                }

                conditionChecker(real4* pos,int* modelId):
                                 pos(pos),modelId(modelId){}

            };

            conditionChecker getConditionChecker(){

                auto pos     = pd->getPos(access::location::gpu, access::mode::read);
                auto model   = pd->getModelId(access::location::gpu, access::mode::read);

                return conditionChecker(pos.raw(),model.raw());
            }
    };

}}}

REGISTER_VERLET_CONDITIONAL_LIST_SET(
    interDifferentType,
    uammd::structured::VerletConditionalListSet<uammd::structured::conditions::interDifferentType>
)
