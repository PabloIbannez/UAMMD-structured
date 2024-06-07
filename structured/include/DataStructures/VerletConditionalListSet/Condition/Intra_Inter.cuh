#ifndef __CONDITION_INTRA_INTER__
#define __CONDITION_INTRA_INTER__

namespace uammd{
namespace structured{
namespace conditions{

    class intra_inter{

        private:

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;

        public:

            /////////////////////////////

            static const int condNum = 2;
            enum condIndex {INTRA=0,INTER=1};

            /////////////////////////////

            intra_inter(std::shared_ptr<GlobalData>            gd,
                        std::shared_ptr<ExtendedParticleData>  pd,
                        DataEntry& dataEntry):gd(gd),pd(pd){

                            System::log<System::MESSAGE>("[Condition] Condition \"intra_inter\" initialized");

                        }

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="inter") {return INTER;}
                if     (conditionName=="intra") {return INTRA;}

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

                int*    modelId;

                inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){}

                inline __device__ void check(const int& i,const int& j,bool cond[condNum]){

                    for(int c=0;c<condNum;c++){
                        cond[c] = false;
                    }

                    if(modelId[i] != modelId[j]){
                        cond[INTER]=true;
                    } else {
                        cond[INTRA]=true;
                    }
                }

                conditionChecker(int* modelId):
                                 modelId(modelId){}
            };

            conditionChecker getConditionChecker(){

                auto model = pd->getModelId(access::location::gpu, access::mode::read);

                return conditionChecker(model.raw());
            }



    };


}}}

#endif
