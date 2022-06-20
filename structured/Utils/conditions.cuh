#ifndef __CONDITIONS__
#define __CONDITIONS__

namespace uammd{
namespace structured{
namespace conditions{
    
    template<class Topology>
    struct all{

        static const int condNum = 1;
        enum condIndex {ALL=0};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="all")     {return ALL;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;

        struct Parameters {};

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;
            return param;
        }

        all(std::shared_ptr<ParticleData>  pd,
            std::shared_ptr<Topology>     top,
            Parameters param):
            sys(pd->getSystem()),pd(pd),
            top(top)
           {
             sys->log<System::MESSAGE>("[Condition] Condition \"all\" initialized");
           }

        
        all(std::shared_ptr<ParticleData>  pd,
            std::shared_ptr<Topology>     top,
            InputFile&                     in):all(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){}

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                cond[0]=true;
            }

            conditionChecker(){}

        };

        conditionChecker getConditionChecker(){
            return conditionChecker();
        }

        size_t getSharedSize(){
            return 0;
        }

    };
    
    template<class Topology>
    struct allInter{

        static const int condNum = 1;
        enum condIndex {INTER=0};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="inter")     {return INTER;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;

        struct Parameters {};

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;
            return param;
        }

        allInter(std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<Topology>     top,
                 Parameters param):
                 sys(pd->getSystem()),pd(pd),
                 top(top)
                {
                  sys->log<System::MESSAGE>("[Condition] Condition \"allInter\" initialized");
                }

        
        allInter(std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<Topology>     top,
                 InputFile&                     in):allInter(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){}

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                cond[0]=true;
            }

            conditionChecker(){}

        };

        conditionChecker getConditionChecker(){
            
            return conditionChecker();
        }

        size_t getSharedSize(){
            return 0;
        }

    };
    
    template<class Topology>
    struct excluded{

        static const int condNum = 1;
        enum condIndex {NONEXCL=0};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="nonExcluded") {return NONEXCL;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excluded(std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<Topology>     top,
                 Parameters param):
                 sys(pd->getSystem()),pd(pd),
                 top(top),
                 exclusionsLabel(param.exclusionsLabel)
                {
                  sys->log<System::MESSAGE>("[Condition] Condition \"excluded\" initialized");

                  exclusionList = std::make_shared<exclusions>(this->pd);
                  exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                }

        
        excluded(std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<Topology>     top,
                 InputFile&                     in):excluded(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            int*  id;

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
                    cond[NONEXCL] = true;
                }
            }

            conditionChecker(int* id,
                             particleExclusionList prtExclList,
                             int maxExclusions):
                             id(id),
                             prtExclList(prtExclList),
                             maxExclusions(maxExclusions){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto id      = pd->getId(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };

    template<class Topology>
    struct chargedExcluded{

        static const int condNum = 2;
        enum condIndex {CHRG=0,NONEXCL=1};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="charged")     {return CHRG;}
            if     (conditionName=="nonExcluded") {return NONEXCL;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        chargedExcluded(std::shared_ptr<ParticleData>  pd,
                        std::shared_ptr<Topology>     top,
                        Parameters param):
                        sys(pd->getSystem()),pd(pd),
                        top(top),
                        exclusionsLabel(param.exclusionsLabel)
                       {
                         sys->log<System::MESSAGE>("[Condition] Condition \"chargedExcluded\" initialized");

                         exclusionList = std::make_shared<exclusions>(this->pd);
                         exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                       }

        
        chargedExcluded(std::shared_ptr<ParticleData>  pd,
                        std::shared_ptr<Topology>     top,
                        InputFile&                     in):chargedExcluded(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            int*  id;
            int*  modelId;
            real* charge;

            particleExclusionList prtExclList;
            int maxExclusions;

            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
            }

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                
                for(int c=0;c<condNum;c++){
                    cond[c] = false;
                }
                
                if(modelId[i] != modelId[j]){return;}
                
                if(abs(charge[i]*charge[j]) > std::numeric_limits<real>::min()){
                    cond[CHRG] = true;
                }
                    
                if(!prtExclList.isPartExcluded(id[j])){
                    cond[NONEXCL] = true;
                }
            }

            conditionChecker(int* id,
                             int* modelId,
                             real* charge,
                             particleExclusionList prtExclList,
                             int maxExclusions):
                             id(id),modelId(modelId),
                             charge(charge),
                             prtExclList(prtExclList),
                             maxExclusions(maxExclusions){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto id      = pd->getId(access::location::gpu, access::mode::read);
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            auto charge  = pd->getCharge(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),
                                    modelId.raw(),
                                    charge.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct excludedIntra{

        static const int condNum = 1;
        enum condIndex {INTRA=0};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")  {return INTRA;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;

        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedIntra(std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<Topology>     top,
                      Parameters param):
                      sys(pd->getSystem()),pd(pd),
                      top(top),
                      exclusionsLabel(param.exclusionsLabel)
                     {
                       sys->log<System::MESSAGE>("[Condition] Condition \"excludedIntra\" initialized");

                       exclusionList = std::make_shared<exclusions>(this->pd);
                       exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                     }

        
        excludedIntra(std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<Topology>     top,
                      InputFile&                     in):excludedIntra(pd,top,inputFileToParam(in)){}
        
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

                if(modelId[i] == modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[INTRA]=true;
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
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct excludedInter{

        static const int condNum = 1;
        enum condIndex {INTER=0};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="inter")  {return INTER;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedInter(std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<Topology>     top,
                      Parameters param):
                      sys(pd->getSystem()),pd(pd),
                      top(top),
                      exclusionsLabel(param.exclusionsLabel)
                     {
                       sys->log<System::MESSAGE>("[Condition] Condition \"excludedInter\" initialized");

                       exclusionList = std::make_shared<exclusions>(this->pd);
                       exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                     }

        
        excludedInter(std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<Topology>     top,
                      InputFile&                     in):excludedInter(pd,top,inputFileToParam(in)){}

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

                if(modelId[i] != modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
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
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct intraInter{

        static const int condNum = 2;
        enum condIndex {INTRA=0,INTER=1};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")  {return INTRA;}
            if     (conditionName=="inter")  {return INTER;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {};

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;
            return param;
        }
        
        intraInter(std::shared_ptr<ParticleData>  pd,
                   std::shared_ptr<Topology>     top,
                   Parameters param):
                   sys(pd->getSystem()),pd(pd),
                   top(top)
                  {
                    sys->log<System::MESSAGE>("[Condition] Condition \"intraInter\" initialized");
                  }

        
        intraInter(std::shared_ptr<ParticleData>  pd,
                   std::shared_ptr<Topology>     top,
                   InputFile&                     in):intraInter(pd,top,inputFileToParam(in)){}


        struct conditionChecker{
            
            int*  modelId;

            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){}

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                
                for(int c=0;c<condNum;c++){
                    cond[c] = false;
                }

                if(modelId[i] == modelId[j]){
                    cond[INTRA]=true;
                } else {
                    cond[INTER]=true;
                }
            }

            conditionChecker(int* modelId):modelId(modelId){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            return conditionChecker(modelId.raw());
        }

        size_t getSharedSize(){
            return 0;
        }

    };
    
    template<class Topology>
    struct excludedIntraInter{

        static const int condNum = 2;
        enum condIndex {INTRA=0,INTER=1};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")  {return INTRA;}
            if     (conditionName=="inter")  {return INTER;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedIntraInter(std::shared_ptr<ParticleData>  pd,
                           std::shared_ptr<Topology>     top,
                           Parameters param):
                           sys(pd->getSystem()),pd(pd),
                           top(top),
                           exclusionsLabel(param.exclusionsLabel)
                          {
                            sys->log<System::MESSAGE>("[Condition] Condition \"excludedIntraInter\" initialized");

                            exclusionList = std::make_shared<exclusions>(this->pd);
                            exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                          }

        
        excludedIntraInter(std::shared_ptr<ParticleData>  pd,
                           std::shared_ptr<Topology>     top,
                           InputFile&                     in):excludedIntraInter(pd,top,inputFileToParam(in)){}

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

                if(modelId[i] == modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[INTRA]=true;
                    } 
                } else {
                    cond[INTER]=true;
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
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct excludedIntraInterCharged{

        static const int condNum = 3;
        enum condIndex {INTRA=0,INTER=1,CHARGED=2};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")  {return INTRA;}
            if     (conditionName=="inter")  {return INTER;}
            if     (conditionName=="charged"){return CHARGED;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedIntraInterCharged(std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<Topology>     top,
                                  Parameters param):
                                  sys(pd->getSystem()),pd(pd),
                                  top(top),
                                  exclusionsLabel(param.exclusionsLabel)
                                 {
                                   sys->log<System::MESSAGE>("[Condition] Condition \"excludedIntraInterCharged\" initialized");

                                   exclusionList = std::make_shared<exclusions>(this->pd);
                                   exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                                 }

        
        excludedIntraInterCharged(std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<Topology>     top,
                                  InputFile&                     in):excludedIntraInterCharged(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            int*  id;
            int*  modelId;

            real* charge;

            particleExclusionList prtExclList;
            int maxExclusions;

            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
            }

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                
                for(int c=0;c<condNum;c++){
                    cond[c] = false;
                }
                
                if(charge[i]*charge[j] != real(0.0)){
                    cond[CHARGED]=true;
                }

                if(modelId[i] == modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[INTRA]=true;
                    } 
                } else {
                    cond[INTER]=true;
                }
            }

            conditionChecker(int* id,int* modelId,real* charge,
                             particleExclusionList prtExclList,
                             int maxExclusions):
                             id(id),modelId(modelId),charge(charge),
                             prtExclList(prtExclList),
                             maxExclusions(maxExclusions){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto id      = pd->getId(access::location::gpu, access::mode::read);
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            auto charge = pd->getCharge(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),charge.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct excludedIntraInterSASAThreshold{

        static const int condNum = 2;
        enum condIndex {INTRA=0,INTER=1};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")  {return INTRA;}
            if     (conditionName=="inter")  {return INTER;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
            
        real SASAThreshold;
        
        struct Parameters {
            real SASAThreshold;
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;
            
            in.getOption("SASAThreshold",InputFile::Required) >> param.SASAThreshold;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedIntraInterSASAThreshold(std::shared_ptr<ParticleData>  pd,
                                        std::shared_ptr<Topology>     top,
                                        Parameters param):
                                        sys(pd->getSystem()),pd(pd),
                                        top(top),
                                        exclusionsLabel(param.exclusionsLabel)
                                       {
                                         sys->log<System::MESSAGE>("[Condition] Condition \"excludedIntraInterSASAThreshold\" initialized");

                                         exclusionList = std::make_shared<exclusions>(this->pd);
                                         exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                                       }

        
        excludedIntraInterSASAThreshold(std::shared_ptr<ParticleData>  pd,
                                        std::shared_ptr<Topology>     top,
                                        InputFile&                     in):excludedIntraInterSASAThreshold(pd,top,inputFileToParam(in)){}

        void setSASAThreshold(real newSASAThreshold){
            SASAThreshold = newSASAThreshold;
        }

        struct conditionChecker{
            
            int*  id;
            int*  modelId;
            
            real* SASA;
            
            real SASAThreshold;

            particleExclusionList prtExclList;
            int maxExclusions;

            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
            }

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                
                for(int c=0;c<condNum;c++){
                    cond[c] = false;
                }

                if(modelId[i] == modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[INTRA]=true;
                    } 
                } else {
                    if(SASA[i] > SASAThreshold and SASA[j] > SASAThreshold){
                        cond[INTER]=true;
                    }
                }
            }

            conditionChecker(int* id,int* modelId,
                             real* SASA,
                             real SASAThreshold,
                             particleExclusionList prtExclList,
                             int maxExclusions):
                             id(id),modelId(modelId),
                             SASA(SASA),
                             SASAThreshold(SASAThreshold),
                             prtExclList(prtExclList),
                             maxExclusions(maxExclusions){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto id      = pd->getId(access::location::gpu, access::mode::read);
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);

            auto SASA    = pd->getSASA(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),
                                    SASA.raw(),
                                    SASAThreshold,
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };
    
    template<class Topology>
    struct excludedIntraInterChargedInter{

        static const int condNum = 3;
        enum condIndex {INTRA=0,INTER=1,CHARGED=2};
        
        int getConditionIndexOf(std::string& conditionName){

            if     (conditionName=="intra")       {return INTRA;}
            if     (conditionName=="inter")       {return INTER;}
            if     (conditionName=="chargedInter"){return CHARGED;}
            else   {
                sys->log<System::CRITICAL>("[Condition]  Requested a condition that is not present, %s",conditionName.c_str());
                return -1;
            }

        }
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData> pd;
        std::shared_ptr<Topology>     top;
        
        struct Parameters {
            std::string exclusionsLabel = "EXCLUSIONS";
        };

        static Parameters inputFileToParam(InputFile& in){
            Parameters param;

            if(in.getOption("EXCLUSIONS_LABEL",InputFile::Optional)){
                in.getOption("EXCLUSIONS_LABEL",InputFile::Optional) >> param.exclusionsLabel;
            }
            return param;
        }
        
        std::string exclusionsLabel;
        std::shared_ptr<exclusions>   exclusionList;

        excludedIntraInterChargedInter(std::shared_ptr<ParticleData>  pd,
                                       std::shared_ptr<Topology>     top,
                                       Parameters param):
                                       sys(pd->getSystem()),pd(pd),
                                       top(top),
                                       exclusionsLabel(param.exclusionsLabel)
                                      {
                                        sys->log<System::MESSAGE>("[Condition] Condition \"excludedIntraInterChargedInter\" initialized");

                                        exclusionList = std::make_shared<exclusions>(this->pd);
                                        exclusionList->loadExclusionListFromTopology(this->top,"EXCLUSIONS");

                                      }

        
        excludedIntraInterChargedInter(std::shared_ptr<ParticleData>  pd,
                                       std::shared_ptr<Topology>     top,
                                       InputFile&                     in):excludedIntraInterChargedInter(pd,top,inputFileToParam(in)){}

        struct conditionChecker{
            
            int*  id;
            int*  modelId;

            real* charge;

            particleExclusionList prtExclList;
            int maxExclusions;

            inline __device__ void set(const int& i,const int& offsetBufferIndex, const char* sharedBuffer){
                prtExclList.set(id[i],(int*)sharedBuffer+offsetBufferIndex*maxExclusions);
            }

            inline __device__ void check(const int& i,const int& j,bool cond[condNum]){
                
                for(int c=0;c<condNum;c++){
                    cond[c] = false;
                }
                
                if(modelId[i] == modelId[j]){
                    if(!prtExclList.isPartExcluded(id[j])){
                        cond[INTRA]=true;
                    } 
                } else {
                    cond[INTER]=true;

                    if(charge[i]*charge[j] != real(0.0)){
                        cond[CHARGED]=true;
                    }

                }
            }

            conditionChecker(int* id,int* modelId,real* charge,
                             particleExclusionList prtExclList,
                             int maxExclusions):
                             id(id),modelId(modelId),charge(charge),
                             prtExclList(prtExclList),
                             maxExclusions(maxExclusions){}
        
        };

        conditionChecker getConditionChecker(){
            
            auto id      = pd->getId(access::location::gpu, access::mode::read);
            auto modelId = pd->getModelId(access::location::gpu, access::mode::read);
            
            auto charge = pd->getCharge(access::location::gpu, access::mode::read);
            
            return conditionChecker(id.raw(),modelId.raw(),charge.raw(),
                                    exclusionList->getParticleExclusionList(),
                                    exclusionList->getMaxExclusions());
        }

        size_t getSharedSize(){
            return exclusionList->getSharedSize();
        }

    };

}}}

#endif
