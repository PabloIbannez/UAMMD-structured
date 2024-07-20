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

    class nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup : public excludedConditionBase {

            std::vector<int> typeGroup1_h;
            std::vector<int> typeGroup2_h;

            thrust::device_vector<int> typeGroup1_d;
            thrust::device_vector<int> typeGroup2_d;

        public:

            /////////////////////////////

            static const int condNum = 4;
            enum condIndex {TYPE_GROUP1_INTRA=0,
                            TYPE_GROUP2_INTRA=1,
                            INTER=2,
                            NOGROUP=3};

            /////////////////////////////

            nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup(std::shared_ptr<GlobalData>            gd,
                                                                                      std::shared_ptr<ExtendedParticleData>  pd,
                                                                                      DataEntry& data):excludedConditionBase(gd,pd,data){

                std::vector<std::string> typeGroup1 = data.getParameter<std::vector<std::string>>("typeGroup1");
                std::vector<std::string> typeGroup2 = data.getParameter<std::vector<std::string>>("typeGroup2");

                auto typeParamHandler = gd->getTypes();

                for(auto& type : typeGroup1){
                    typeGroup1_h.push_back(typeParamHandler->getTypeId(type));
                }

                for(auto& type : typeGroup2){
                    typeGroup2_h.push_back(typeParamHandler->getTypeId(type));
                }

                // Sort typeGroup1 and typeGroup2
                std::sort(typeGroup1_h.begin(),typeGroup1_h.end());
                std::sort(typeGroup2_h.begin(),typeGroup2_h.end());

                typeGroup1_d = typeGroup1_h;
                typeGroup2_d = typeGroup2_h;

                System::log<System::MESSAGE>("[Condition] Condition "
                                             "\"nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExcInter_nonExclNoGroup\" "
                                             "initialized.");

            }

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="typeGroup1Intra"){return TYPE_GROUP1_INTRA;}
                else if(conditionName=="typeGroup2Intra"){return TYPE_GROUP2_INTRA;}
                else if(conditionName=="inter")          {return INTER;}
                else if(conditionName=="noGroup")        {return NOGROUP;}
                else
                {
                    System::log<System::CRITICAL>("[Condition] Requested a condition"
                                                  " that is not present, %s",
                                                  conditionName.c_str());
                    return -1;
                }

                return -1;

            }

            ////////////////////////////////////////////

            struct conditionChecker{

                int* id;
                int* modelId;

                real4*  pos;

                int* typeGroup1;
                int* typeGroup2;

                int typeGroup1Size;
                int typeGroup2Size;

                particleExclusionList prtExclList;
                int maxExclusions;

                __device__ int binarySearch(const int* data, int size, int value) {
                    int left = 0;
                    int right = size - 1;

                    while (left <= right) {
                        int middle = left + (right - left) / 2;

                        if (data[middle] == value) {
                            return middle; // or return a specific flag indicating found
                        }
                        if (data[middle] < value) {
                            left = middle + 1;
                        }
                        else {
                            right = middle - 1;
                        }
                    }

                    return -1; // or return a specific flag indicating not found
                }

                inline __device__ int getGroup(int type){

                    if     (binarySearch(typeGroup1,typeGroup1Size,type)!=-1){
                        return 1;
                    }
                    else if(binarySearch(typeGroup2,typeGroup2Size,type)!=-1){
                        return 2;
                    }
                    else{
                        return 0;
                    }

                    return 0;

                }

                inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                    prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
                }

                inline __device__ void check(const int& i,const int& j,bool cond[condNum]){

                    for(int c=0;c<condNum;c++){
                        cond[c] = false;
                    }

                    if(!prtExclList.isPartExcluded(id[j])){

                        int modelId_i = modelId[i];
                        int modelId_j = modelId[j];

                        int type_i = int(pos[i].w);
                        int type_j = int(pos[j].w);

                        int group_i = getGroup(type_i);
                        int group_j = getGroup(type_j);

                        bool noGroup = (group_i==0 || group_j==0);

                        if(noGroup){
                            cond[NOGROUP] = true;
                        }else{
                            if (modelId_i != modelId_j){
                                cond[INTER] = true;
                            }
                            else{
                                if(group_i!=group_j){
                                    cond[INTER] = true;
                                } else{
                                    if(group_i==1){
                                        cond[TYPE_GROUP1_INTRA] = true;
                                    }
                                    else{
                                        cond[TYPE_GROUP2_INTRA] = true;
                                    }
                                }
                            }
                        }
                    }
                }

                conditionChecker(int* id, int* modelId,
                                 real4* pos,
                                 int* typeGroup1,int* typeGroup2,
                                 int typeGroup1Size,int typeGroup2Size,
                                 particleExclusionList prtExclList,
                                 int maxExclusions):
                                 id(id),modelId(modelId),
                                 pos(pos),
                                 typeGroup1(typeGroup1),typeGroup2(typeGroup2),
                                 typeGroup1Size(typeGroup1Size),typeGroup2Size(typeGroup2Size),
                                 prtExclList(prtExclList),
                                 maxExclusions(maxExclusions){}

            };


            conditionChecker getConditionChecker(){

                int*      id = pd->getId(access::location::gpu,access::mode::read).raw();
                int* modelId = pd->getModelId(access::location::gpu,access::mode::read).raw();

                real4* pos = pd->getPos(access::location::gpu,access::mode::read).raw();

                int* typeGroup1 = thrust::raw_pointer_cast(typeGroup1_d.data());
                int* typeGroup2 = thrust::raw_pointer_cast(typeGroup2_d.data());

                return conditionChecker(id,modelId,
                                        pos,
                                        typeGroup1,typeGroup2,
                                        typeGroup1_d.size(),typeGroup2_d.size(),
                                        exclusionList->getParticleExclusionList(),
                                        exclusionList->getMaxExclusions());
            }
    };

}}}

REGISTER_VERLET_CONDITIONAL_LIST_SET(
    nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup,
    uammd::structured::VerletConditionalListSet<uammd::structured::conditions::nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup>
)
