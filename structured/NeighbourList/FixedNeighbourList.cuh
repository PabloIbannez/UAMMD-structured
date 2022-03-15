#ifndef __FIXED_NEIGHBOUR_LIST__
#define __FIXED_NEIGHBOUR_LIST__

namespace uammd{
namespace structured{
  
    namespace FixedNeighbourList_ns{
    
    __global__ void id2indexLists(const int* neighbourList_id,
                                        int* neighbourList_index,
                                  const int* neighbourStart,
                                  const int* nNeighbourPerParticle,
                                  const int* id2index,
                                  int nEntries){
            
            int index = blockIdx.x*blockDim.x + threadIdx.x;
            if(index>=nEntries) return;

            const int* neigBegin       = neighbourList_id+index;
                  int* neigBegin_index = neighbourList_index+index;

            for(int n=0;n<nNeighbourPerParticle[index];n++){
                neigBegin_index[n*nEntries]=id2index[neigBegin[n*nEntries]];
            }
        }
    
    }
    
    template<class Topology>
    class FixedNeighbourList : public ParameterUpdatable{

        private:
            
            std::shared_ptr<System>      sys;
            std::shared_ptr<ParticleData> pd;
            std::shared_ptr<ParticleGroup> pg;
            std::shared_ptr<Topology>    top;           
                
            std::string topologyLabel;

            int nEntries;
            int maxNeighbours;

            thrust::host_vector<int> h_neighbourList;
            thrust::host_vector<int> h_neighbourStart;
            thrust::host_vector<int> h_nNeighbourPerParticle;
            
            thrust::device_vector<int> neighbourList_id;
            thrust::device_vector<int> neighbourList_index;
            
            thrust::device_vector<int> neighbourStart;
            thrust::device_vector<int> nNeighbourPerParticle;
         
            bool needsUpdate = true;

            connection reorderConnection;
          
            void handleReorder(){
                sys->log<System::DEBUG1>("[FixedNeighboutList] Issuing a list update after a reorder.");
                needsUpdate = true;
            }

        public:
            
            struct Parameters{
                std::string topologyLabel;
            };
            
            struct NeighbourListData{
                int N;
                const int * neighbourList;
                const int * neighbourStart;
                const int * numberNeighbours;
            };

            FixedNeighbourList(std::shared_ptr<System>       sys,
                               std::shared_ptr<ParticleData>  pd,
                               std::shared_ptr<ParticleGroup> pg,
                               std::shared_ptr<Topology>     top,
                               Parameters param):sys(sys),
                                                 pd(pd),pg(pg),
                                                 top(top),
                                                 topologyLabel(param.topologyLabel){
              
                reorderConnection  = pd->getReorderSignal()->connect([this](){this->handleReorder();});

                sys->log<System::MESSAGE>("[FixedNeighbourList] Initialized");
                
                std::map<int,std::vector<int>> neigL;
                top->loadNeighbourList(topologyLabel,neigL);
                
                nEntries      = 0;
                maxNeighbours = 0;

                for(auto const& [id,neigPerPartL] : neigL){
                    nEntries ++;
                    int cNeig = neigPerPartL.size();
                    if(cNeig > maxNeighbours){
                        maxNeighbours = cNeig;
                    }
                }
                
                sys->log<System::MESSAGE>("[FixedNeighbourList] Number of entries: %i, "
                                           "max number of neighbours: %i",
                                           nEntries,maxNeighbours);
                
                //Load data to inner representation
                h_neighbourList.resize(nEntries*(maxNeighbours+1),-1);  //+1 for self reference at [0]
                h_neighbourStart.resize(nEntries);
                h_nNeighbourPerParticle.resize(nEntries,0);
            
                neighbourList_index.resize(nEntries*(maxNeighbours+1)); //+1 for self reference at [0]
                
                uint index=0;
                for(auto const& [id,neigPerPartL] : neigL){
                    h_neighbourStart[index]=index;
                    int neigCounter=0;
                    h_neighbourList[neigCounter*nEntries+index]=id;
                    neigCounter++;
                    for(int n : neigPerPartL){
                        h_neighbourList[neigCounter*nEntries+index]=n;
                        neigCounter++;
                    }
                    h_nNeighbourPerParticle[index]=neigCounter;

                    index++;
                }
                    
                neighbourList_id      = h_neighbourList;
                nNeighbourPerParticle = h_nNeighbourPerParticle;
                neighbourStart        = h_neighbourStart;
            }
          
            ~FixedNeighbourList(){
                reorderConnection.disconnect();
            }

            int getMaxNeighbours(){
                return maxNeighbours;
            }

            void update(cudaStream_t st = 0){
                if(needsUpdate){
                    
                    auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                    const int * neigList_id = thrust::raw_pointer_cast(neighbourList_id.data());
                    int * neigList_index    = thrust::raw_pointer_cast(neighbourList_index.data());
                    
                    const int * neigStart  = thrust::raw_pointer_cast(neighbourStart.data());
                    const int * nNeigPPart = thrust::raw_pointer_cast(nNeighbourPerParticle.data());
                    
                    int Nthreads=128;
                    int Nblocks=nEntries/Nthreads + ((nEntries%Nthreads)?1:0);

                    FixedNeighbourList_ns::id2indexLists<<<Nthreads,Nblocks,0,st>>>(neigList_id,
                                                                                    neigList_index,
                                                                                    neigStart,
                                                                                    nNeigPPart,
                                                                                    id2index,
                                                                                    nEntries);
                    
                    CudaCheckError();

                    needsUpdate=false;
                }

            }
            
            NeighbourListData getNeighbourList(std::string conditionName,cudaStream_t st = 0){
                
                NeighbourListData nl;
                    
                nl.N                = nEntries; 
                nl.neighbourList    = thrust::raw_pointer_cast(neighbourList_index.data());
                nl.numberNeighbours = thrust::raw_pointer_cast(nNeighbourPerParticle.data());
                nl.neighbourStart   = thrust::raw_pointer_cast(neighbourStart.data());
                
                return nl;
            }

    };
}}

#endif
