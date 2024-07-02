#pragma once

#include"utils/Grid.cuh"
#include<thrust/device_vector.h>

#include"Interactor/NeighbourList/CellList/CellListBase.cuh"
#include"Interactor/NeighbourList/CellList/NeighbourContainer.cuh"
#include"Interactor/NeighbourList/BasicList/BasicListBase.cuh"

#include"DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"

namespace uammd{
namespace structured{

    namespace VerletConditionalListSet_ns{

        template<class basicListInfo,class condition>
        __global__ void fillNeighbourConditionalListSet(CellList_ns::NeighbourContainer ni,
                                                        ParticleGroup::IndexIterator group2GlobalIndex,
                                                        basicListInfo* bListInfo,
                                                        int maxNeighboursPerParticle,
                                                        real cutOff2,
                                                        int numberParticles, Box box,
                                                        uint* tooManyNeighboursFlag,
                                                        int* batchId,
                                                        typename condition::conditionChecker condChecker){

            int localIndex = blockIdx.x*blockDim.x + threadIdx.x;

            extern __shared__ char sharedBuffer[];

            if(localIndex>=numberParticles) return;

            const int i_global = group2GlobalIndex[ni.getGroupIndexes()[localIndex]];

            const auto sortPos = ni.getSortedPositions();
            const real3 pi     = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + localIndex));

            ni.set(localIndex);
            condChecker.set(i_global,threadIdx.x,sharedBuffer);

            //Initialization
            int nneigh[condition::condNum];
            for(int c=0;c<condition::condNum;c++){
                nneigh[c]=1;
                bListInfo[c].neighbourList_ptr[localIndex] = i_global;
            }

            bool condArray[condition::condNum];

            auto it = ni.begin();
            while(it){
                auto n = *it++;

                const int j_global = group2GlobalIndex[n.getGroupIndex()];

                if (batchId[i_global]==batchId[j_global]){

                    const real3 pj  = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + n.getInternalIndex()));

                    const real3 rij = box.apply_pbc(pj-pi);
                    const real  r2 = dot(rij, rij);

                    if(r2 <=cutOff2){

                        condChecker.check(i_global,
                                          j_global,condArray);

                        for(int c=0;c<condition::condNum;c++){

                            if(condArray[c]){

                                if(nneigh[c]>=maxNeighboursPerParticle){
                                    nneigh[c]++;
                                    atomicMax(tooManyNeighboursFlag, nneigh[c]);
                                    return;
                                }

                                if(i_global != j_global){ //Avoid self
                                    bListInfo[c].neighbourList_ptr[nneigh[c]*numberParticles + localIndex] = j_global;
                                    nneigh[c]++;
                                }
                            }
                        }
                    }

                }
            }

            for(int c=0;c<condition::condNum;c++){
                bListInfo[c].numberNeighbours_ptr[localIndex]=nneigh[c];
            }

        }

        template<class PosIterator>
        __global__ void checkMaximumDrift(PosIterator currentPos,
            				                      real4 *storedPos,
            				                      real maxDistAllowed,
            				                      uint* errorFlag,
            				                      Box box,
            				                      int numberParticles){

          int id = blockIdx.x*blockDim.x + threadIdx.x;
          if(id>=numberParticles) return;

          real3 currentPos_i  = make_real3(currentPos[id]);
          real3 previousPos_i = make_real3(storedPos[id]);

          real3 rij = box.apply_pbc(currentPos_i - previousPos_i);

          if(dot(rij, rij)>=(maxDistAllowed*maxDistAllowed)){atomicAdd(errorFlag, 1u);}
        }
    }

    template<class condition>
    class VerletConditionalListSet: public VerletConditionalListSetBase {

        public:

            struct Parameters{
              real cutOff = 0.0;
              real cutOffVerletFactor = 1.1;
          };

        private:

            std::string name;

            //Box
            Box box;

            //Cut-off radii

            real cutOff;

            real cutOffVerletFactor;
            real cutOffVerlet;

            real thresholdDistance;

            //Cell list

            CellListBase cl;

            //Lists data

            int maxNeighboursPerParticle=0;
            bool forceNextRebuild=true;

        public:

            //basicListInfo_ptr is public to avoid nvcc complaints, but should not be used directly
            struct basicListInfo_ptr{
                int* neighbourList_ptr;
                int* numberNeighbours_ptr;
            };

        private:

            struct basicListInfo{
                thrust::device_vector<int> neighbourList;
                thrust::device_vector<int> numberNeighbours;
            };

            std::vector<basicListInfo>               neigCondListSet;
            thrust::device_vector<basicListInfo_ptr> neigCondListSet_ptr;

            //Condtion

            std::shared_ptr<condition> cond;

            //Flags

            uint*                        tooManyNeighboursFlagCPU;
            thrust::device_vector<uint>  tooManyNeighboursFlagGPU;

            uint*                        particleDriftFlagCPU;
            thrust::device_vector<uint>  particleDriftFlagGPU;

            //Connection

            connection reorderConnection;

            //Drift

            thrust::device_vector<real4> storedPosition;

            ////////////////////////////////////////////////////////////

            void increaseMaximumNeighboursPerParticle(){

                int N = pg->getNumberParticles();

                this->maxNeighboursPerParticle += 1;

                System::log<System::DEBUG2>("[VerletConditionalListSet] (%s) Increasing maximum number of neighbours to %d",name.c_str(),maxNeighboursPerParticle);

                for(basicListInfo& list : neigCondListSet){
                    list.neighbourList.resize(N*maxNeighboursPerParticle);
                    list.numberNeighbours.resize(N);
                }

                for(int ncl=0;ncl<condition::condNum;ncl++){
                    neigCondListSet_ptr[ncl]={thrust::raw_pointer_cast(neigCondListSet[ncl].neighbourList.data()),
                                              thrust::raw_pointer_cast(neigCondListSet[ncl].numberNeighbours.data())};
                }

            }

            void fillNeighbourConditionalListSet(cudaStream_t st){

                bool filled = false;
                while(not filled){

                    System::log<System::DEBUG3>("[VerletConditionalListSet] (%s) Attempting to fill list with %d neighbours per particle",name.c_str(), maxNeighboursPerParticle);

                    int   N = pg->getNumberParticles();
                    auto  group2GlobalIndex = pg->getIndexIterator(access::location::gpu);

                    auto batchId = pd->getBatchId(access::location::gpu, access::mode::write);

                    auto cldata = cl.getCellList();
                    auto ni = CellList_ns::NeighbourContainer(cldata);

                    CudaCheckError();

                    tooManyNeighboursFlagCPU[0] = 0;
                    cudaMemcpyAsync(thrust::raw_pointer_cast(tooManyNeighboursFlagGPU.data()), tooManyNeighboursFlagCPU,
                                    sizeof(uint), cudaMemcpyHostToDevice,st);

                    uint* d_tooManyNeighboursFlag = thrust::raw_pointer_cast(tooManyNeighboursFlagGPU.data());

                    int Nthreads=256;
                    int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                    VerletConditionalListSet_ns::fillNeighbourConditionalListSet<basicListInfo_ptr,condition>
                    <<<Nblocks, Nthreads,  Nthreads*cond->getSharedSize(), st>>>(ni,
                                                                                 group2GlobalIndex,
                                                                                 thrust::raw_pointer_cast(neigCondListSet_ptr.data()),
                                                                                 maxNeighboursPerParticle,
                                                                                 cutOffVerlet*cutOffVerlet,
                                                                                 N, box,
                                                                                 d_tooManyNeighboursFlag,
                                                                                 batchId.raw(),
                                                                                 cond->getConditionChecker());
                    cudaMemcpyAsync(tooManyNeighboursFlagCPU,
                                    thrust::raw_pointer_cast(tooManyNeighboursFlagGPU.data()),
                                    sizeof(uint), cudaMemcpyDeviceToHost,st);

                    cudaStreamSynchronize(st);
                    uint flag = tooManyNeighboursFlagCPU[0];

                    CudaCheckError();

                    filled = not (flag != 0);

                    if(not filled){
                        increaseMaximumNeighboursPerParticle();
                    }
                }

            }

            void storeCurrentPosition(cudaStream_t st){

                auto pos     = pd->getPos(access::location::gpu, access::mode::read);

                int numberParticles = pg->getNumberParticles();

                storedPosition.resize(numberParticles);

                thrust::copy(thrust::cuda::par.on(st), pos, pos + numberParticles, storedPosition.begin());
                cudaStreamSynchronize(st);
            }

            void handleReorder(){
                System::log<System::DEBUG1>("[ConditionedVerletListSet] Issuing a list update after a reorder.");
                forceNextRebuild = true;
            }

        public:

            VerletConditionalListSet(std::shared_ptr<GlobalData>            gd,
                                     std::shared_ptr<ParticleGroup>         pg,
                                     DataEntry& data,
                                     std::string name):VerletConditionalListSetBase(gd,pg),
                                                       name(name){

                this->cutOffVerletFactor = data.getParameter<real>("cutOffVerletFactor",1.1);
                this->setCutOff(data.getParameter<real>("cutOff",0.0));

                ////////////////////////////

                this->box = gd->getEnsemble()->getBox();

                ////////////////////////////
                //Set up condition

                std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());
                cond = std::make_shared<condition>(gd,pd,data);

                ////////////////////////////

                cudaError_t status = cudaMallocHost((void**)&tooManyNeighboursFlagCPU, sizeof(uint));
                if (status != cudaSuccess){
                    System::log<System::CRITICAL>("[VerletConditionalListSet] (%s) Error allocating pinned host memory",name.c_str());
                }

                tooManyNeighboursFlagGPU.resize(1);

                status = cudaMallocHost((void**)&particleDriftFlagCPU, sizeof(uint));
                if (status != cudaSuccess){
                    System::log<System::CRITICAL>("[VerletConditionalListSet] (%s) Error allocating pinned host memory",name.c_str());
                }

                particleDriftFlagGPU.resize(1);

                ////////////////////////////

                reorderConnection  = pd->getReorderSignal()->connect([this](){this->handleReorder();});

                ////////////////////////////

                neigCondListSet.resize(condition::condNum);
                neigCondListSet_ptr.resize(condition::condNum);

                ////////////////////////////

                this->increaseMaximumNeighboursPerParticle();

            }

            ~VerletConditionalListSet(){
                reorderConnection.disconnect();
                cudaFreeHost(tooManyNeighboursFlagCPU);
                cudaFreeHost(particleDriftFlagCPU);
            }

            std::string getName() override {
                return name;
            }

            void update(cudaStream_t st) override{

                if(not forceNextRebuild){

                    int N = pg->getNumberParticles();

                    auto pos     = pd->getPos(access::location::gpu, access::mode::read);
                    auto posIter = pg->getPropertyIterator(pos);


                    particleDriftFlagCPU[0] = 0;
                    cudaMemcpyAsync(thrust::raw_pointer_cast(particleDriftFlagGPU.data()), particleDriftFlagCPU,
                                    sizeof(uint), cudaMemcpyHostToDevice,st);

                    int Nthreads=128;
                    int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);

                    VerletConditionalListSet_ns::checkMaximumDrift
                    <<<Nblocks, Nthreads, 0, st>>>(posIter,
                                                   thrust::raw_pointer_cast(storedPosition.data()),
                                                   thresholdDistance,
                                                   thrust::raw_pointer_cast(particleDriftFlagGPU.data()),
                                                   box,
                                                   N);
                    cudaMemcpyAsync(particleDriftFlagCPU,
                                    thrust::raw_pointer_cast(particleDriftFlagGPU.data()),
                                    sizeof(uint), cudaMemcpyDeviceToHost,st);
                    CudaCheckError();
                    cudaStreamSynchronize(st);
                    uint particleDriftFlag = particleDriftFlagCPU[0];
                    if(particleDriftFlag>0){
                        System::log<System::DEBUG2>("[VerletList] Found %d particles over threshold.", particleDriftFlag);
                    } else {
                        return;
                    }

                }

                if(forceNextRebuild){
                    forceNextRebuild = false;
                }

                //Perform update

                {
                    //Store current pos
                    storeCurrentPosition(st);

                    //
                    int N = pg->getNumberParticles();

                    auto pos = pd->getPos(access::location::gpu, access::mode::read);
                    auto posIter = pg->getPropertyIterator(pos);

                    //Update cellList
                    Grid grid(box, cutOffVerlet);
                    cl.update(posIter, N, grid, st);

                    fillNeighbourConditionalListSet(st);
                    CudaCheckError();

                }

            }

            NeighbourListData getNeighbourList(std::string conditionName) override{

                ////Write neig list to file
                //std::ofstream out("list.dat");
                //for(int i=0;i<condition::condNum;i++){
                //    thrust::host_vector<int> tmpList = neigCondListSet[i].neighbourList;
                //    thrust::host_vector<int> tmpN    = neigCondListSet[i].numberNeighbours;

                //    out << "list: " << i << std::endl;
                //    for(int j=0;j<pg->getNumberParticles();j++){
                //        out << tmpList[j] << ": ";
                //        for(int n=1;n<tmpN[j];n++){
                //            out << tmpList[n*pg->getNumberParticles()+j] << " ";
                //        }
                //        out << std::endl;
                //    }

                //}
                //out.close();

                //std::cout << "PAUSE" << std::endl;
                //std::cin.get();

                int condIndex =  cond->getConditionIndexOf(conditionName);

                NeighbourListData nl;


                nl.N = pg->getNumberParticles();

                nl.neighbourList       = thrust::raw_pointer_cast(neigCondListSet[condIndex].neighbourList.data());
                nl.numberNeighbours    = thrust::raw_pointer_cast(neigCondListSet[condIndex].numberNeighbours.data());
                nl.neighbourStart      = StrideIterator(0);

                return nl;
            }

            void setCutOff(real newCutOff) override {

                if(cutOff == newCutOff){
                    System::log<System::DEBUG5>("[VerletConditionalListSet] (%s) Cut-off unchanged: %f",name.c_str(),cutOff);
                    return;
                }

                cutOff         = newCutOff;
                cutOffVerlet   = cutOff*cutOffVerletFactor;

                thresholdDistance = (cutOffVerlet-cutOff)/2.0;

                forceNextRebuild = true;

                System::log<System::MESSAGE>("[VerletConditionalListSet] (%s) Cell list cut-off updated: %f, cut-off Verlet : %f",
                                             name.c_str(),
                                             this->cutOff,this->cutOffVerlet);

            }

            real getCutOff() override { return cutOff;}
    };

}}
