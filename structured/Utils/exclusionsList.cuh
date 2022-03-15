#ifndef __EXCLUSION_LIST__
#define __EXCLUSION_LIST__

namespace uammd{
namespace structured{
            
    struct particleExclusionList{

        int** exclusionStart;
        int*  nExclusionPerParticle;
    
        int current_i;
        int current_nExc; 
    
        int* start;
        int* startBuffer;
    
        particleExclusionList(int** exclusionStart,
                              int*  nExclusionPerParticle):exclusionStart(exclusionStart),
                                                           nExclusionPerParticle(nExclusionPerParticle){}
    
        inline __device__ void set(int i){
            current_i=i;
            current_nExc=nExclusionPerParticle[current_i];
            start = exclusionStart[current_i];
    
            startBuffer=start;               
        }
    
        inline __device__ void set(int i,int* strB){
            current_i=i;
            current_nExc=nExclusionPerParticle[current_i];
            start = exclusionStart[current_i];
            startBuffer = strB;
    
            forj(0,current_nExc){
                startBuffer[j]=start[j];               
            }
        }
    
        inline __device__ bool isPartExcluded(int k){
    
            forj(0,current_nExc){
                if(startBuffer[j]==k){return true;}
            }
    
            return false;
    
        }
    
    };

    class exclusions{
        
        private:
            
            std::shared_ptr<System>      sys;
            std::shared_ptr<ParticleData> pd;
                
            int nExclusions;
            int maxExclusions;

            thrust::host_vector<int>  h_exclusionList;
            thrust::host_vector<int*> h_exclusionStart;
            thrust::host_vector<int>  h_nExclusionPerParticle;
            
            thrust::device_vector<int>  exclusionList;
            thrust::device_vector<int*> exclusionStart;
            thrust::device_vector<int>  nExclusionPerParticle;
        
        public:
            
            exclusions(std::shared_ptr<System>      sys,
                       std::shared_ptr<ParticleData> pd):sys(sys),
                                                         pd(pd){}

            template<class Topology>
            void loadExclusionListFromTopology(std::shared_ptr<Topology> top, std::string topologyLabel){

                sys->log<System::MESSAGE>("[Exclusions] Initialized");
                
                int numParticles = pd->getNumParticles();
            
                std::map<int,std::vector<int>> exclL;
                top->loadNeighbourList(topologyLabel,exclL);
                
                nExclusions   = 0;
                maxExclusions = 0;

                for(auto const& [id,exclPerPartL] : exclL){
                    int cExcl = exclPerPartL.size();
                    nExclusions += cExcl;
                    if(cExcl > maxExclusions){
                        maxExclusions = cExcl;
                    }
                }
                
                sys->log<System::MESSAGE>("[Exclusions] Number of exclusions: %i, "
                                           "max number of exclusions: %i",
                                           nExclusions,maxExclusions);
                
                //Load data to inner representation
                h_exclusionList.resize(nExclusions);
                h_exclusionStart.resize(numParticles,nullptr);
                h_nExclusionPerParticle.resize(numParticles,0);
                
                int exclCounter=0;
                for(auto const& [id,exclPerPartL] : exclL){
                    for(int e : exclPerPartL){
                        h_exclusionList[exclCounter]=e;
                        exclCounter++;
                    }
                    h_nExclusionPerParticle[id]=exclPerPartL.size();
                }
                    
                exclusionList         = h_exclusionList;
                nExclusionPerParticle = h_nExclusionPerParticle;
                    
                int offSet=0;
                fori(0,numParticles){
                    if(h_nExclusionPerParticle[i]!=0){
                        h_exclusionStart[i]=thrust::raw_pointer_cast(exclusionList.data())+offSet;
                        offSet=offSet+h_nExclusionPerParticle[i];
                    }
                }
                    
                exclusionStart = h_exclusionStart;
            }

            int getMaxExclusions(){
                return maxExclusions;
            }
            
            size_t getSharedSize(){
                return this->getMaxExclusions()*sizeof(int);
            }

            particleExclusionList getParticleExclusionList(){
                return particleExclusionList(thrust::raw_pointer_cast(exclusionStart.data()),
                                             thrust::raw_pointer_cast(nExclusionPerParticle.data()));
            }
    };
}}

#endif
