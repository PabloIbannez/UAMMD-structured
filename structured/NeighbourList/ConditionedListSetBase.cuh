#ifndef CONDITIONEDLISTSETBASE_CUH
#define CONDITIONEDLISTSETBASE_CUH

#include"utils/Grid.cuh"
#include<thrust/device_vector.h>
#include"Interactor/NeighbourList/CellList/CellListBase.cuh"
#include"Interactor/NeighbourList/CellList/NeighbourContainer.cuh"
#include"Interactor/NeighbourList/BasicList/BasicListBase.cuh"

namespace uammd{
namespace structured{
    
    namespace ConditionedNeighbourListSet_ns{
        
        template<class basicListInfo,class condition>
        __global__ void fillConditionedNeighbourListSet(CellList_ns::NeighbourContainer ni,
                                                        ParticleGroup::IndexIterator group2GlobalIndex,
                                                        basicListInfo* bListInfo,
                                                        int maxNeighboursPerParticle,
                                                        real cutOff2,
                                                        int numberParticles, Box box,
                                                        uint* tooManyNeighboursFlag,
                                                        int* simId,
                                                        typename condition::conditionChecker condChecker){
            int id = blockIdx.x*blockDim.x + threadIdx.x;

            extern __shared__ char sharedBuffer[];

            if(id>=numberParticles) return;

            const int i_group  = ni.getGroupIndexes()[id];
            const int i_global = group2GlobalIndex[i_group];

            const auto sortPos = ni.getSortedPositions();
            const real3 pi     = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + id));

            ni.set(id);
            condChecker.set(i_global,threadIdx.x,sharedBuffer);

            //Initialization
            int nneigh[condition::condNum];
            for(int c=0;c<condition::condNum;c++){
                nneigh[c]=1;
                bListInfo[c].neighbourList_ptr[id] = i_global;
            }

            bool condArray[condition::condNum];

            auto it = ni.begin();
            while(it){
                auto n = *it++;
                const int j_group  = n.getGroupIndex();
                const int j_global = group2GlobalIndex[j_group];
                const real3 pj  = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + n.getInternalIndex()));
                const real3 rij = box.apply_pbc(pj-pi);
                const real r2 = dot(rij, rij);
                if(r2 <=cutOff2 and 
                   simId[i_global]==simId[j_global]){
                    condChecker.check(i_global,
                                      j_global,condArray);
                    for(int c=0;c<condition::condNum;c++){
                        if(condArray[c]){
                            if(nneigh[c]>=maxNeighboursPerParticle){
                                nneigh[c]++;
                                atomicMax(tooManyNeighboursFlag, nneigh[c]);
                                return;
                            }
                            if(i_group != j_group){ //Avoid self
                                bListInfo[c].neighbourList_ptr[nneigh[c]*numberParticles + id] = j_global;
                                nneigh[c]++;
                            }
                        }
                    }
                }
            }

            for(int c=0;c<condition::condNum;c++){
                bListInfo[c].numberNeighbours_ptr[id]=nneigh[c];
            }

        }
    }

    template<class condition>
    class ConditionedListSetBase{
        
        protected:
            
            /*//Row-major
            struct NeighbourListOffsetFunctor{
                NeighbourListOffsetFunctor(int str):
                    stride(str){}
                int stride;
                inline __host__ __device__ int operator()(const int &index) const{
                    return index*stride;
                }
            };*/
          
            shared_ptr<ParticleData> pd;
            shared_ptr<ParticleGroup> pg;
            shared_ptr<System> sys;

            CellListBase cl;

            struct basicListInfo{
                thrust::device_vector<int> neighbourList; 
                thrust::device_vector<int> numberNeighbours;
            };

            uint*                        errorFlagsCPU;
            thrust::device_vector<uint>  errorFlagsGPU;

            int maxNeighboursPerParticle=0;

            std::shared_ptr<condition> cond;

            /*//Row-major
            using CountingIterator = cub::CountingInputIterator<int>;
            using StrideIterator   = cub::TransformInputIterator<int, NeighbourListOffsetFunctor, CountingIterator>;
            */
            using StrideIterator = cub::CountingInputIterator<int>;
        
        public:
          
            struct basicListInfo_ptr{
                int* neighbourList_ptr; 
                int* numberNeighbours_ptr;
            };
    
            struct NeighbourListData{
                int   N;
                const int * neighbourList;
                const int * numberNeighbours;
                const int*  internal2GroupIndex;
                /*//Row-major
                StrideIterator neighbourStart = StrideIterator(CountingIterator(0), NeighbourListOffsetFunctor(0));
                */
                StrideIterator neighbourStart = StrideIterator(0);
            };

            std::vector<basicListInfo>               condNeigListSet;
            thrust::device_vector<basicListInfo_ptr> condNeigListSet_ptr;

            ConditionedListSetBase(shared_ptr<System> sys,
                                   shared_ptr<ParticleData> pd,
                                   shared_ptr<ParticleGroup> pg,
                                   shared_ptr<condition> cond):sys(sys), pd(pd), pg(pg),
                                                cond(cond){
                cudaError_t status = cudaMallocHost((void**)&errorFlagsCPU, sizeof(uint));
                if (status != cudaSuccess){
                    sys->log<System::CRITICAL>("[ConditionedListSetBase] Error allocating pinned host memory");
                }
                
                errorFlagsGPU.resize(1);

                condNeigListSet.resize(condition::condNum);
                condNeigListSet_ptr.resize(condition::condNum);

                this->increaseMaximumNeighboursPerParticle();
            }

            ~ConditionedListSetBase(){
                cudaFreeHost(errorFlagsCPU);
            }

            template<class PositionIterator>
            void update(PositionIterator pos, int numberParticles, Box box, real cutOff, cudaStream_t st = 0){
                
                int N = pg->getNumberParticles();
                Grid grid(box, cutOff);
                cl.update(pos, N, grid, st);
                fillConditionedNeighbourListSet(box,cutOff,st);
                CudaCheckError();
            }
            
            NeighbourListData getNeighbourList(std::string conditionName,cudaStream_t st = 0){

                int condIndex =  cond->getConditionIndexOf(conditionName);
                
                NeighbourListData nl;
                
                auto cld = cl.getCellList();

                nl.N = pg->getNumberParticles();
                nl.neighbourList = thrust::raw_pointer_cast(condNeigListSet[condIndex].neighbourList.data());
                nl.numberNeighbours = thrust::raw_pointer_cast(condNeigListSet[condIndex].numberNeighbours.data());
                /*//Row-major
                nl.neighbourStart = StrideIterator(CountingIterator(0), NeighbourListOffsetFunctor(maxNeighboursPerParticle));
                */
                nl.internal2GroupIndex = cld.groupIndex;
                nl.neighbourStart = StrideIterator(0);
                return nl;
            }

        private:
            
            void increaseMaximumNeighboursPerParticle(){
                
                int N = pg->getNumberParticles();
                
                this->maxNeighboursPerParticle += 32;
                
                System::log<System::MESSAGE>("[ConditionedNeighbourListSet] Increasing maximum number of neighbours to %d", maxNeighboursPerParticle);
                
                for(basicListInfo& list : condNeigListSet){
                    list.neighbourList.resize(N*maxNeighboursPerParticle);
                    list.numberNeighbours.resize(N);
                }
                
                for(int ncl=0;ncl<condition::condNum;ncl++){
                    condNeigListSet_ptr[ncl]={thrust::raw_pointer_cast(condNeigListSet[ncl].neighbourList.data()),
                                              thrust::raw_pointer_cast(condNeigListSet[ncl].numberNeighbours.data())};
                }
                
            }

            void fillConditionedNeighbourListSet(Box box,real cutOff,cudaStream_t st){
                while(not tryToFillConditionedNeighbourListSet(box,cutOff,st)){
                    increaseMaximumNeighboursPerParticle();
                }
            }

            bool tryToFillConditionedNeighbourListSet(Box box,real cutOff,cudaStream_t st){
                System::log<System::DEBUG3>("[ConditionedNeighbourListSet] Attempting to fill list with %d neighbours per particle", maxNeighboursPerParticle);
                
                int   N = pg->getNumberParticles();
                auto  group2GlobalIndex = pg->getIndexIterator(access::location::gpu);

                auto simId = pd->getSimulationId(access::location::gpu, access::mode::write);
                
                auto cldata = cl.getCellList();
                auto ni = CellList_ns::NeighbourContainer(cldata);
                
                errorFlagsCPU[0] = 0;
                cudaMemcpyAsync(thrust::raw_pointer_cast(errorFlagsGPU.data()), errorFlagsCPU, 
                                sizeof(uint), cudaMemcpyHostToDevice,st);

                uint* d_tooManyNeighboursFlag = thrust::raw_pointer_cast(errorFlagsGPU.data());

                int Nthreads=128;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                ConditionedNeighbourListSet_ns::fillConditionedNeighbourListSet<basicListInfo_ptr,condition>
                <<<Nblocks, Nthreads,  Nthreads*cond->getSharedSize(), st>>>(ni,
                                                                             group2GlobalIndex,
                                                                             thrust::raw_pointer_cast(condNeigListSet_ptr.data()),
                                                                             maxNeighboursPerParticle,
                                                                             cutOff*cutOff,
                                                                             N, box,
                                                                             d_tooManyNeighboursFlag,
                                                                             simId.raw(),
                                                                             cond->getConditionChecker());
                cudaMemcpyAsync(errorFlagsCPU,
                                thrust::raw_pointer_cast(errorFlagsGPU.data()),  
                                sizeof(uint), cudaMemcpyDeviceToHost,st);

                cudaStreamSynchronize(st);
                uint flag = errorFlagsCPU[0];

                bool foundTooManyNeighbours = flag != 0;
                return not foundTooManyNeighbours;
            }
    };

}}
#endif

