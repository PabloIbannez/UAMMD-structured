#ifndef __SELECTORS__
#define __SELECTORS__

namespace uammd{
namespace structured{
namespace selectors{
    
    class idMax{

        private:
            
            int maxId;

        public:
            
            idMax(int maxId):maxId(maxId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partId = (pd->getId(access::cpu, access::read).raw())[particleIndex];

                return (partId<maxId);
            }


    };
    
    class idSet{

        private:
            
            std::vector<int> ids;

        public:
            
            idSet(std::vector<int>& ids):ids(ids){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partId = (pd->getId(access::cpu, access::read).raw())[particleIndex];

                for(int& id : ids){
                    if(partId == id){
                        return true;
                    }
                }

                return false;
            }


    };
    
    class id{

        private:
            
            int i;

        public:
            
            id(int i):i(i){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partId = (pd->getId(access::cpu, access::read).raw())[particleIndex];
                return (partId==i);
            }


    };
    
    class notId{

        private:
            
            int i;

        public:
            
            notId(int i):i(i){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partId = (pd->getId(access::cpu, access::read).raw())[particleIndex];
                return (partId!=i);
            }


    };
    
    class chainId{

        private:
            
            int chId;

        public:
            
            chainId(int chId):chId(chId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partChId = (pd->getChainId(access::cpu, access::read).raw())[particleIndex];
                return (partChId==chId);
            }


    };
    
    class modelId{

        private:
            
            int mdlId;

        public:
            
            modelId(int mdlId):mdlId(mdlId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partMdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];
                return (partMdlId==mdlId);
            }
    };

    class modelIdChainId{

        private:
            
            int mdlId;
            int chnId;

        public:
            
            modelIdChainId(int mdlId,
                           int chnId):mdlId(mdlId),
                                      chnId(chnId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partMdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];
                int partChnId = (pd->getChainId(access::cpu, access::read).raw())[particleIndex];
                return (partMdlId==mdlId and partChnId==chnId);
            }
    };
    
    class notModelId{

        private:
            
            int mdlId;

        public:
            
            notModelId(int mdlId):mdlId(mdlId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partMdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];
                return (partMdlId!=mdlId);
            }


    };
    
    class simulationId{

        private:
            
            int simId;

        public:
            
            simulationId(int simId):simId(simId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partSimId = (pd->getSimulationId(access::cpu, access::read).raw())[particleIndex];
                return (partSimId==simId);
            }


    };
    
    class simulationIdModelId{

        private:
            
            int simId;
            int mdlId;

        public:
            
            simulationIdModelId(int simId,
                                int mdlId):simId(simId),
                                           mdlId(mdlId){};

            bool isSelected(int particleIndex, std::shared_ptr<ParticleData> &pd){
                int partSimId = (pd->getSimulationId(access::cpu, access::read).raw())[particleIndex];
                int partMdlId = (pd->getModelId(access::cpu, access::read).raw())[particleIndex];
                return (partSimId==simId and partMdlId==mdlId);
            }


    };
    

}}}

#endif
