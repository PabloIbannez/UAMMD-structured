#ifndef __CONDITION_ALL__
#define __CONDITION_ALL__

namespace uammd{
namespace structured{
namespace conditions{

    class all{

        private:

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;

        public:

            /////////////////////////////

            static const int condNum = 1;
            enum condIndex {ALL=0};

            /////////////////////////////

            all(std::shared_ptr<GlobalData>            gd,
                std::shared_ptr<ExtendedParticleData>  pd,
                DataEntry& dataEntry):gd(gd),pd(pd){

                    System::log<System::MESSAGE>("[Condition] Condition \"all\" initialized");

                }

            ////////////////////////////////////////////

            int getConditionIndexOf(std::string& conditionName){

                if     (conditionName=="all") {return ALL;}

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

                inline __device__ void set(const int& i,
                                           const int& offsetBufferIndex,
                                           const char* sharedBuffer){}

                inline __device__ void check(const int& i,const int& j,
                                             bool cond[condNum]){
                    cond[0]=true;
                }

                conditionChecker(){}
            };


            conditionChecker getConditionChecker(){
                return conditionChecker();
            }

    };

}}}

#endif
