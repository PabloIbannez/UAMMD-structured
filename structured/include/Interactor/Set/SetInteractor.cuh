#pragma once

#include"Interactor/Interactor.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace SetInteractor_ns{

        template<class transverser>
        struct reduction{

            using type = typename transverser::reduceType;

            __device__ __forceinline__ type operator () (const type& reduced1,
                                                         const type& reduced2){

                return transverser::reduce(reduced1,reduced2);
            }
        };

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser>
        __device__ __forceinline__ void performThreadReduction(const int* id2index,
                                                               const int& groupIndex,
                                                               const int& N,const int* start,
                                                               ComputationalData computational,
                                                               SetParameters*  __restrict__ setParam,
                                                               PotentialTransverser& transverser,
                                                               const int& reductionStep,
                                                               const int& setIndex,
                                                               typename PotentialTransverser::reduceType reduced[nReductions][nSet]){


            typename PotentialTransverser::reduceType redAcc = transverser.zero(computational,setParam[groupIndex],reductionStep,setIndex,reduced);

            for(int i=0;i<N;i++){
                const int index = id2index[start[i]];
                redAcc = transverser.reduce(redAcc,transverser.transform(index,computational,setParam[groupIndex],reductionStep,setIndex,reduced));
            }

            reduced[reductionStep][setIndex]=redAcc;

        }

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser>
        __device__ __forceinline__ void spreadAndSetThread(const int* id2index,
                                                           const int& groupIndex,
                                                           const int& N,const int* start,
                                                           ComputationalData computational,
                                                           SetParameters*  __restrict__ setParam,
                                                           PotentialTransverser& transverser,
                                                           typename PotentialTransverser::reduceType reduced[nReductions][nSet],
                                                           const typename PotentialTransverser::propertyType& property,
                                                           const int& setIndex){


            for(int i=0;i<N;i++){
                const int index = id2index[start[i]];
                transverser.set(index,
                                transverser.spread(property,
                                                   index,
                                                   computational,
                                                   setParam[groupIndex],
                                                   setIndex,
                                                   reduced),
                                computational,
                                setParam[groupIndex],
                                setIndex,
                                reduced);
            }
        }

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser>
        __global__ void transverseSet_ThreadPerSet(const int*  __restrict__ id2index,
                                                         int** __restrict__ listStart[nSet],
                                                         int*  __restrict__ n[nSet],
                                                         int nGroups,
                                                   ComputationalData computational,
                                                   SetParameters*  __restrict__ setParam,
                                                   PotentialTransverser transverser) {

            const int groupIndex = blockIdx.x*blockDim.x + threadIdx.x;

            if(groupIndex>=nGroups){return;}

            typename PotentialTransverser::reduceType reduced[nReductions][nSet];

            for(int redStep = 0; redStep < nReductions; redStep++){
                for(int setIndex = 0; setIndex < nSet; setIndex++){

                    const int* start = listStart[setIndex][groupIndex];
                    const int  N     = n[setIndex][groupIndex];

                    performThreadReduction
                    <nReductions,nSet,
                    ComputationalData,SetParameters,
                    PotentialTransverser>
                    (id2index,
                     groupIndex,
                     N,start,
                     computational,
                     setParam,
                     transverser,
                     redStep,
                     setIndex,
                     reduced);
                }
            }

            typename PotentialTransverser::propertyType property = transverser.compute(computational,setParam[groupIndex],reduced);

            for(int setIndex = 0; setIndex < nSet; setIndex++){

                const int* start = listStart[setIndex][groupIndex];
                const int  N     = n[setIndex][groupIndex];

                spreadAndSetThread
                <nReductions,nSet,
                 ComputationalData,SetParameters,
                 PotentialTransverser>
                (id2index,
                 groupIndex,
                 N,start,
                 computational,
                 setParam,
                 transverser,
                 reduced,
                 property,
                 setIndex);
            }

        }

        //BLOCK

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser, int THREADS_PER_BLOCK,
                 typename BlockReduce>
        __device__ __forceinline__ void performBlockReduction(const int* id2index,
                                                              const int& groupIndex,const int& particleIndex,
                                                              const int& N,const int* start,
                                                              ComputationalData computational,
                                                              SetParameters*  __restrict__ setParam,
                                                              PotentialTransverser& transverser,
                                                              typename BlockReduce::TempStorage* temp_storage,
                                                              const int& reductionStep,
                                                              const int& setIndex,
                                                              typename PotentialTransverser::reduceType reduced[nReductions][nSet]){


                typename PotentialTransverser::reduceType redAcc = transverser.zero(computational,
                                                                                    setParam[groupIndex],
                                                                                    reductionStep,setIndex,
                                                                                    reduced);
                typename PotentialTransverser::reduceType redBuffer;

                const int nMax = N/THREADS_PER_BLOCK + ((N%THREADS_PER_BLOCK)?1:0);
                for(int n=0;n<nMax;n++){

                    if(n>0){
                        __syncthreads(); //All threads enter here
                    }

                    int maxThread;
                    if((n+1)*blockDim.x<N){
                        maxThread = blockDim.x;
                    } else {
                        maxThread = N-n*blockDim.x;
                    }

                    int i=particleIndex+n*blockDim.x;

                    typename PotentialTransverser::reduceType data = transverser.zero(computational,setParam[groupIndex],
                                                                                      reductionStep,
                                                                                      setIndex,reduced);

                    if(i<N){
                        const int index = id2index[start[i]];
                        data = transverser.transform(index,
                                                     computational,setParam[groupIndex],
                                                     reductionStep,setIndex,reduced);
                    }

                    redBuffer = BlockReduce(*temp_storage).Reduce(data,reduction<PotentialTransverser>(),maxThread);

                    if(threadIdx.x==0){
                        redAcc+=redBuffer;
                    }
                }

                __syncthreads(); //Not needed really?
                if(threadIdx.x==0){
                    reduced[reductionStep][setIndex]=redAcc;
                }
                __syncthreads();
        }

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser>
        __device__ __forceinline__ void spreadAndSetBlock(const int* id2index,
                                                          const int& groupIndex,const int& particleIndex,
                                                          const int& N,const int* start,
                                                          ComputationalData computational,
                                                          SetParameters*  __restrict__ setParam,
                                                          PotentialTransverser& transverser,
                                                          typename PotentialTransverser::reduceType reduced[nReductions][nSet],
                                                          const typename PotentialTransverser::propertyType& property,
                                                          const int& setIndex){


            for(int i=particleIndex;i<N;i+=blockDim.x){
                const int index = id2index[start[i]];
                transverser.set(index,
                                transverser.spread(property,
                                                   index,
                                                   computational,
                                                   setParam[groupIndex],
                                                   setIndex,
                                                   reduced),
                                computational,
                                setParam[groupIndex],
                                setIndex,
                                reduced);
            }
        }

        template<int nReductions, int nSet,
                 class ComputationalData, class SetParameters,
                 class PotentialTransverser, int THREADS_PER_BLOCK>
        __global__ void transverseSet_BlockPerSet(const int*  __restrict__ id2index,
                                                        int** __restrict__ listStart[nSet],
                                                        int*  __restrict__ n[nSet],
                                                  ComputationalData computational,
                                                  SetParameters*  __restrict__ setParam,
                                                  PotentialTransverser transverser){

            //
            typedef cub::BlockReduce<typename PotentialTransverser::reduceType, THREADS_PER_BLOCK> BlockReduce;

            __shared__ typename BlockReduce::TempStorage temp_storage;
            //

            //
            __shared__ typename PotentialTransverser::reduceType reduced[nReductions][nSet];

            const int groupIndex    = blockIdx.x;
            const int particleIndex = threadIdx.x;

            //Block reduction
            for(int redStep = 0; redStep < nReductions; redStep++){
                for(int setIndex = 0; setIndex < nSet; setIndex++){

                    const int* start = listStart[setIndex][groupIndex];
                    const int  N     = n[setIndex][groupIndex]; //Particles in group

                    performBlockReduction
                    <nReductions,nSet,
                     ComputationalData,SetParameters,
                     PotentialTransverser,THREADS_PER_BLOCK,BlockReduce>
                    (id2index,
                     groupIndex,particleIndex,
                     N,start,
                     computational,
                     setParam,
                     transverser,
                     &temp_storage,
                     redStep,
                     setIndex,
                     reduced);
                }
            }

            typename PotentialTransverser::propertyType property = transverser.compute(computational,setParam[groupIndex],reduced);

            for(int setIndex = 0; setIndex < nSet; setIndex++){

                const int* start = listStart[setIndex][groupIndex];
                const int  N     = n[setIndex][groupIndex]; //Particles in group

                spreadAndSetBlock
                <nReductions,nSet,
                 ComputationalData,SetParameters,
                 PotentialTransverser>
                (id2index,
                 groupIndex,particleIndex,
                 N,start,
                 computational,
                 setParam,
                 transverser,
                 reduced,
                 property,setIndex);
            }

        }

    }

    template<class PotentialType, int SET_SIZE_THRESHOLD = 32, int THREADS_PER_BLOCK=512>
    class SetInteractor: public Interactor {

        private:

            using SetParameters = typename PotentialType::SetParameters;

            ////////////////////////////////////

            int nSet = PotentialType::nSet;

            std::shared_ptr<GlobalData>     gd;
            std::shared_ptr<PotentialType>  pot;

            ////////////////////////////////////

            bool large = false;

            ////////////////////////////////////

            std::vector<std::shared_ptr<groupsList>> sets;
            thrust::device_vector<SetParameters>     setParam_d;

            thrust::device_vector<int**> listStart_d;
            thrust::device_vector<int*>  listSize_d;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy        = false;
            bool warningForce         = false;

            void sumSet_ThreadPerSet(Computables comp,cudaStream_t st){

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<PotentialType>::value){

                        int*** listStart = thrust::raw_pointer_cast(listStart_d.data());
                        int**  n         = thrust::raw_pointer_cast(listSize_d.data());

                        auto setParam_ptr  = thrust::raw_pointer_cast(setParam_d.data());

                        int  nGroups  = setParam_d.size();

                        const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads = THREADS_PER_BLOCK;
                        int Nblocks  = nGroups/Nthreads + ((nGroups%Nthreads)?1:0);

                        SetInteractor_ns::transverseSet_ThreadPerSet
                        <PotentialType::nReductions,
                         PotentialType::nSet,
                         typename PotentialType::ComputationalData,
                         typename PotentialType::SetParameters,
                         typename PotentialType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(id2index,
                                                     listStart,
                                                     n,
                                                     nGroups,
                                                     pot->getComputationalData(comp,st),
                                                     setParam_ptr,
                                                     pot->getEnergyTransverser());

                        CudaCheckError();
                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[SetInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }
                }

                if(comp.force == true){

                    if constexpr (has_getForceTransverser<PotentialType>::value){

                        int*** listStart = thrust::raw_pointer_cast(listStart_d.data());
                        int**  n         = thrust::raw_pointer_cast(listSize_d.data());

                        auto setParam_ptr  = thrust::raw_pointer_cast(setParam_d.data());

                        int  nGroups  = setParam_d.size();

                        const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads = THREADS_PER_BLOCK;
                        int Nblocks  = nGroups/Nthreads + ((nGroups%Nthreads)?1:0);

                        SetInteractor_ns::transverseSet_ThreadPerSet
                        <PotentialType::nReductions,
                         PotentialType::nSet,
                         typename PotentialType::ComputationalData,
                         typename PotentialType::SetParameters,
                         typename PotentialType::ForceTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(id2index,
                                                     listStart,
                                                     n,
                                                     nGroups,
                                                     pot->getComputationalData(comp,st),
                                                     setParam_ptr,
                                                     pot->getForceTransverser());

                        CudaCheckError();
                    } else {
                        if(!warningForce){
                            System::log<System::WARNING>("[SetInteractor] (%s) Requested non-implemented transverser (force)",
                                                         name.c_str());
                            warningForce = true;
                        }
                    }
                }

            }

            void sumSet_BlockPerSet(Computables comp,cudaStream_t st){

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<PotentialType>::value){

                        int*** listStart = thrust::raw_pointer_cast(listStart_d.data());
                        int**  n         = thrust::raw_pointer_cast(listSize_d.data());

                        auto setParam_ptr  = thrust::raw_pointer_cast(setParam_d.data());

                        int  nGroups  = setParam_d.size();

                        const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nGroups;

                        SetInteractor_ns::transverseSet_BlockPerSet
                        <PotentialType::nReductions,
                         PotentialType::nSet,
                         typename PotentialType::ComputationalData,
                         typename PotentialType::SetParameters,
                         typename PotentialType::EnergyTransverser,
                         THREADS_PER_BLOCK>
                        <<<Nblocks,Nthreads,0, st>>>(id2index,
                                                     listStart,
                                                     n,
                                                     pot->getComputationalData(comp,st),
                                                     setParam_ptr,
                                                     pot->getEnergyTransverser());
                        CudaCheckError();
                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[SetInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }

                }

                if(comp.force == true){

                    if constexpr (has_getForceTransverser<PotentialType>::value){

                        int*** listStart = thrust::raw_pointer_cast(listStart_d.data());
                        int**  n         = thrust::raw_pointer_cast(listSize_d.data());

                        auto setParam_ptr  = thrust::raw_pointer_cast(setParam_d.data());

                        int  nGroups  = setParam_d.size();

                        const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nGroups;

                        SetInteractor_ns::transverseSet_BlockPerSet
                        <PotentialType::nReductions,
                         PotentialType::nSet,
                         typename PotentialType::ComputationalData,
                         typename PotentialType::SetParameters,
                         typename PotentialType::ForceTransverser,
                         THREADS_PER_BLOCK>
                        <<<Nblocks,Nthreads,0, st>>>(id2index,
                                                     listStart,
                                                     n,
                                                     pot->getComputationalData(comp,st),
                                                     setParam_ptr,
                                                     pot->getForceTransverser());
                        CudaCheckError();
                    } else {
                        if(!warningForce){
                            System::log<System::WARNING>("[SetInteractor] (%s) Requested non-implemented transverser (force)",
                                                         name.c_str());
                            warningForce = true;
                        }
                    }

                }
            }

        public:

            SetInteractor(std::shared_ptr<GlobalData>           gd,
                          std::shared_ptr<ParticleGroup>        pg,
                          DataEntry& data,
                          std::shared_ptr<PotentialType> pot,
                          std::string name):Interactor(pg,"SetInteractor: \"" +name+"\""),
                                            gd(gd),
                                            pot(pot){


                ////////////////////////////

                thrust::host_vector<SetParameters>  setParam_h;

                thrust::host_vector<int*>   listSize_h;
                thrust::host_vector<int**>  listStart_h;

                auto setLabels = pot->getSetLabels();

                //Print the given labels, order matters!!
                {
                    std::string labels = "";
                    for(int i = 0; i < nSet; i++){
                      if(i == nSet - 1){
                        labels += setLabels[i] + " (index: " + std::to_string(i) + ")";
                      }else{
                        labels += setLabels[i] + " (index: " + std::to_string(i) + ") ----";
                      }
                    }
                    System::log<System::MESSAGE>("[SetInteractor] (%s) Sets: %s",name.c_str(),labels.c_str());
                }

                {
                    std::vector<std::vector<std::vector<int>>> sets2check;

                    for(int i=0;i<nSet;i++){

                        System::log<System::MESSAGE>("[SetInteractor] (%s) Loading set \"%s\" ...",name.c_str(),setLabels[i].c_str());

                        //Load set
                        auto set_raw = data.getData<std::vector<int>>(setLabels[i]);
                        //Check set not intersects
                        {
                            sets2check.push_back(set_raw);

                            std::vector<int> intersection = setsIntersection(sets2check);

                            if(intersection.size() != 0){
                                std::string message;
                                for(auto i : intersection){
                                    message+=std::string(" ")+std::to_string(i);
                                }
                                System::log<System::CRITICAL>("[SetInteractor] (%s) Some elements (%s) appear in different ids groups",
                                                               name.c_str(),message.c_str());
                            }
                        }

                        //Load set from set_raw
                        sets.push_back(std::make_shared<groupsList>(set_raw));

                        System::log<System::MESSAGE>("[SetInteractor] (%s) Loaded %d ids sets from \"%s\"",name.c_str(),set_raw.size(),setLabels[i].c_str());
                    }
                }

                {
                    System::log<System::MESSAGE>("[SetInteractor] (%s) Loading set parameters ...",name.c_str());

                    //Load set parameters
                    std::vector<std::string> ignoredLabels = setLabels;
                    auto setParametersMap = data.getDataMapIgnoringLabels(ignoredLabels);

                    for(auto& i : setParametersMap){
                        setParam_h.push_back(pot->processSetParameters(i));
                    }

                    //Copy to gpu
                    setParam_d = setParam_h;

                    //Check groupslist and setParameters have the same size
                    for(int i=0;i<nSet;i++){
                        if(this->sets[i]->getNGroups() != setParam_d.size()){
                            System::log<System::CRITICAL>("[SetInteractor] (%s) The number of given groups (%i) and the parameters given (%i) do not match for set %i",
                                                           name.c_str(),this->sets[i]->getNGroups(),setParam_d.size(),i);
                        }
                    }

                    for(int i=0;i<nSet;i++){
                        if(this->sets[i]->getMaxGroupSize() > SET_SIZE_THRESHOLD){
                            large=true;
                            System::log<System::MESSAGE>("[SetInteractor] (%s) Group %i (which belongs to \"%s\")has got %i particles."
                                                         " It is considered large (>%i).",
                                                          name.c_str(),i,setLabels[i].c_str(),this->sets[i]->getMaxGroupSize(),SET_SIZE_THRESHOLD);
                        }
                    }

                    if(large){
                        System::log<System::MESSAGE>("[SetInteractor] (%s) At least one group is considered large. (>%i).",
                                                      name.c_str(),SET_SIZE_THRESHOLD);
                        System::log<System::MESSAGE>("[SetInteractor] (%s) One block per set algorithm will be used."
                                                     " Threads per block: %i.",name.c_str(),THREADS_PER_BLOCK);
                    } else {
                        System::log<System::MESSAGE>("[SetInteractor] (%s) All groups are considered small. (<= %i)",name.c_str(),SET_SIZE_THRESHOLD);
                        System::log<System::MESSAGE>("[SetInteractor] (%s) One thread per set algorithm will be used.",name.c_str());
                    }

                }

                //Store sets start and size
                {
                    listStart_h.resize(nSet);
                    listSize_h.resize(nSet);

                    for(int i=0;i<nSet;i++){
                        auto setParam  = sets[i]->getGroupsListInfo(access::location::gpu);

                        listStart_h[i] = setParam.listStart;
                        listSize_h[i]  = setParam.n;
                    }

                    listStart_d = listStart_h;
                    listSize_d  = listSize_h;
                }


            }

            std::shared_ptr<PotentialType> getPotential(){
                return pot;
            }

            void sum(Computables comp,cudaStream_t st) override {

                if(large){
                    sumSet_BlockPerSet(comp,st);
                } else {
                    sumSet_ThreadPerSet(comp,st);
                }

            }
    };

}}}
