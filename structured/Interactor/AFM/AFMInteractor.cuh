#ifndef __AFM_INTERACTOR__
#define __AFM_INTERACTOR__

#include"Interactor/Interactor.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace AFMInteractor_ns {

        template<class ComputationalData, class AFMParameters,class AFMTransverser>
        __global__ void transverseAFM(const int nAFM,
                                      const int nParticlesSample,
                                      const int* __restrict__    tipIds,
                                      const int* __restrict__ sampleIds,
                                      const int* __restrict__  id2index,
                                      const int* __restrict__     batId,
                                      const ComputationalData           computational,
                                      const AFMParameters* __restrict__ afmParameters,
                                      const int*  __restrict__ id2BufferIndex,
                                      typename AFMTransverser::resultType*  __restrict__ quantityBuffer,
                                      AFMTransverser afmTransverser,
                                      bool isForce) {

            const int indexLocal = blockIdx.x*blockDim.x + threadIdx.x;

            if(indexLocal==0){
                //Perform operations over the tip
                for(int i=int(0);i<nAFM;i++){
                    const int tipId    = tipIds[i];
                    const int tipIndex = id2index[tipId];

                    typename AFMTransverser::resultType quantityTip = afmTransverser.zero();
                    afmTransverser.accumulate(quantityTip,afmTransverser.compute(tipIndex,
                                                                                 tipIndex,
                                                                                 computational,afmParameters[i]));

                    //Note we not set the value, we add to the buffer
                    const int bufferTipIndex = id2BufferIndex[tipId];
                    quantityBuffer[bufferTipIndex] = quantityTip + afmTransverser.get(tipIndex);
                    //Note we add the storaged value since other interactions may have been applied over the tip
                }
            }

            if(indexLocal>=nParticlesSample){return;}

            {//Sample
                const int id          = sampleIds[indexLocal];
                const int indexGlobal = id2index[id];
                const int bId         = batId[indexGlobal];

                const int tipId       = tipIds[bId];
                const int tipIndex    = id2index[tipId];

                typename AFMTransverser::resultType quantityLocal = afmTransverser.zero();

                afmTransverser.accumulate(quantityLocal,afmTransverser.compute(indexGlobal,
                                                                               tipIndex,
                                                                               computational,afmParameters[bId]));
                afmTransverser.set(indexGlobal,quantityLocal);

                const int bufferIndex = id2BufferIndex[id];

                //Note not +=, we are setting the value. The force is the opposite of the quantity.
                if( isForce ){
                    quantityBuffer[bufferIndex] = -quantityLocal;
                } else {
                    quantityBuffer[bufferIndex] =  quantityLocal;
                }
            }

        }

    }

    template<class AFMType>
    class AFMInteractor: public Interactor{

        private:

            const int THREADS_PER_BLOCK = 256;

            using AFMParameters = typename AFMType::AFMParameters;

            int nAFM;
            int nParticlesSample;

            std::shared_ptr<GlobalData> gd;
            std::shared_ptr<AFMType>   afm;

            ////////////////////////////////////

            std::shared_ptr<groupsList> tip;
            std::shared_ptr<groupsList> sample;

            std::vector<std::string> labels = {"idSet_i", "idSet_j"}; // idSet_i = tip, idSet_j = sample
            //We use this labels to ensure that everythig works fine with pyUAMMD

            thrust::host_vector<AFMParameters>   afmParam_h;
            thrust::device_vector<AFMParameters> afmParam_d;

            ////////////////////////////////////
            //Buffering

            thrust::device_vector<int> id2BufferIndex_d;
            thrust::device_vector<int> bufferOffset_d;

            thrust::device_vector<real>   energyBuffer_d;
            thrust::device_vector<real4>  forceBuffer_d;

            ////////////////////////////////////
            //Reduction

            void*  cubReductionTempStorage    ;

            size_t cubEnergyReductionTempStorageSize;
            size_t cubForceReductionTempStorageSize;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy = false;
            bool warningForce  = false;

        public:

            ////////////////////////////////////

            //We need to create an iterator for reduction output
            //This has to send the data to the correct tip
            //Example:
            //auto tipParam    = tip->getGroupsListInfo(access::location::gpu);
            //Reduction index, tip id, tip index,
            //0              , tipParam.list[0] , id2index[tipParam.list[0]],
            //1              , tipParam.list[1] , id2index[tipParam.list[1]],
            //...

            struct reductionIndex2TipIndex: public thrust::unary_function<int,int>{

                const int* __restrict__ id2index;
                const int* __restrict__ tipIds;

                const int nAFM;

                reductionIndex2TipIndex(const int* __restrict__ tipIds,
                                        const int* __restrict__ id2index,
                                        const int nAFM):
                    tipIds(tipIds),
                    id2index(id2index),
                    nAFM(nAFM){}

                __host__ __device__ int operator()(const int& index) const {
                    if(index<nAFM){
                        return id2index[tipIds[index]];
                    }else{
                        return -1;
                    }
                }
            };


            using reductionIndex2TipIndexIteratorTransform   = thrust::transform_iterator<reductionIndex2TipIndex,
                                                                                          thrust::counting_iterator<int>>;

            using energyPermutationIterator = thrust::permutation_iterator<real*,
                                                                           reductionIndex2TipIndexIteratorTransform>;

            using forcePermutationIterator  = thrust::permutation_iterator<real4*,
                                                                           reductionIndex2TipIndexIteratorTransform>;

            energyPermutationIterator getEnergyPermutationIterator(){

                auto id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                reductionIndex2TipIndex transFunctor(tip->getGroupsListInfo(access::location::gpu).list,
                                                     id2index,
                                                     nAFM);

                reductionIndex2TipIndexIteratorTransform transIterator = thrust::make_transform_iterator(thrust::counting_iterator<int>(0),
                                                                                                         transFunctor);


                auto energy = pd->getEnergy(access::location::gpu,access::mode::readwrite);

                return thrust::make_permutation_iterator(energy.begin(),transIterator);
            }

            forcePermutationIterator getForcePermutationIterator(){

                auto id2index = this->pd->getIdOrderedIndices(access::location::gpu);

                reductionIndex2TipIndex transFunctor(tip->getGroupsListInfo(access::location::gpu).list,
                                                     id2index,
                                                     nAFM);

                reductionIndex2TipIndexIteratorTransform transIterator = thrust::make_transform_iterator(thrust::counting_iterator<int>(0),
                                                                                                         transFunctor);

                auto force = pd->getForce(access::location::gpu,access::mode::readwrite);

                return thrust::make_permutation_iterator(force.begin(),transIterator);
            }


            //Each data row is considered to be an indepenent afm.
            //For each afm we expect to found two sets of particles
            //one for the tip and one for the sample:
            // -The size of the tip is expected to be 1
            // -The size of the sample is not determined, it could be 0 or more.
            //
            //If the simulation is made of more than one batch, we expect to found
            //an afm for each batch. All batches must have one afm defined.
            //The batch id must match the afm id.
            //
            //Example:
            //"labels": ["idSet_i", "idSet_j"], // idSet_i = tip, idSet_j = sample
            //"data": [
            // [0,[1,2,3,4,5,...,100]], // afm 0, the batch id for particles 0,1,2,3,4,5,...,100 must be 0
            // [101,[102,103,104,105,...,250]], // afm 1, the batch id for particles 102,103,104,105,...,250 must be 1
            // ...
            // ]

            AFMInteractor(std::shared_ptr<GlobalData>    gd,
                          std::shared_ptr<ParticleGroup> pg,
                          DataEntry& data,
                          std::shared_ptr<AFMType> afm,
                          std::string name):Interactor(pg,"AFMInteractor: \"" +name+"\""),
                                            gd(gd),
                                            afm(afm){

                ////////////////////////////

                //Print idSet_i = tip, idSet_j = sample
                System::log<System::MESSAGE>("[AFMInteractor] (%s) Tip (%s), Sample (%s)",
                                              name.c_str(),
                                              labels[0].c_str(), labels[1].c_str());

                {
                    std::vector<std::vector<std::vector<int>>> sets2check;

                    // Load tip
                    System::log<System::MESSAGE>("[AFMInteractor] Loading tip ...");
                    std::vector<std::vector<int>> tip_raw = data.getData<std::vector<int>>(labels[0]);
                    //Check all tips have only one particle
                    for(auto& tip_i: tip_raw){
                        if(tip_i.size() == 0){
                            System::log<System::CRITICAL>("[AFMInteractor] (%s) Tip set is empty!",name.c_str());
                        }
                        if(tip_i.size() > 1){
                            System::log<System::CRITICAL>("[AFMInteractor] (%s) Only single particle tips are supported!",name.c_str());
                        }
                    }

                    sets2check.push_back(tip_raw);
                    tip = std::make_shared<groupsList>(tip_raw);

                    // Load sample
                    System::log<System::MESSAGE>("[AFMInteractor] Loading sample ...");
                    std::vector<std::vector<int>> sample_raw = data.getData<std::vector<int>>(labels[1]);
                    sample = std::make_shared<groupsList>(sample_raw);
                    sets2check.push_back(sample_raw);

                    //Check set not intersects
                    {
                        std::vector<int> intersection = SetInteractor_ns::setsIntersection(sets2check);

                        if(intersection.size() != 0){
                            std::string message;
                            for(auto i : intersection){
                                message+=std::string(" ")+std::to_string(i);
                            }
                            System::log<System::CRITICAL>("[AFMInteractor] (%s) Some elements (%s) appear in different ids groups",
                                                           name.c_str(),message.c_str());
                        }
                    }

                    {
                        auto batId = this->pd->getBatchId(access::location::cpu,
                                                               access::mode::read);
                        auto id2index = this->pd->getIdOrderedIndices(access::location::cpu);

                        auto tip_h    = tip->getGroupsListInfo(access::location::cpu);
                        auto sample_h = sample->getGroupsListInfo(access::location::cpu);

                        std::set<int> presentBatchIds;
                        for(int i = 0; i < tip->getNGroups(); ++i){
                            int tipId    = tip_h.listStart[i][0];
                            int tipBatchId = batId[id2index[tipId]];
                            //Check group index (i) is equal to the batch id
                            if(tipBatchId != i){
                                System::log<System::CRITICAL>("[AFMInteractor] (%s) Tip group index (%d) is not equal to the batch id (%d)",
                                                               name.c_str(),i,tipBatchId);
                            }
                            presentBatchIds.insert(tipBatchId);
                            for(int j = 0; j < sample_h.n[i]; ++j){
                                int sampleId    = sample_h.listStart[i][j];
                                int sampleBatchId = batId[id2index[sampleId]];
                                if(tipBatchId != sampleBatchId){
                                    System::log<System::CRITICAL>("[AFMInteractor] (%s) Tip and sample must be in the same batch!."
                                                                  "But particle with id %d is in batch %d and tip %d is in batch %d",
                                                                   name.c_str(),
                                                                   sampleId,sampleBatchId,
                                                                   tipId,tipBatchId);


                                }
                            }
                        }

                        std::set<int> totalBatchIds;
                        for(int i = 0; i < pd->getNumParticles(); ++i){
                            totalBatchIds.insert(batId[i]);
                        }

                        //Check that all batches are present
                        if(presentBatchIds.size() != totalBatchIds.size()){
                            std::string message;
                            for(auto i : totalBatchIds){
                                if(presentBatchIds.find(i) == presentBatchIds.end()){
                                    message+=std::string(" ")+std::to_string(i);
                                }
                            }
                            System::log<System::CRITICAL>("[AFMInteractor] (%s) Some batches (%s) are not present in tip or sample",
                                                           name.c_str(),message.c_str());
                        }

                    }

                    nAFM = tip->getNGroups();
                    nParticlesSample = sample->getNRelations();
                }

                {
                    System::log<System::MESSAGE>("[AFMInteractor] (%s) Loading afm parameters ...",name.c_str());

                    //Load set parameters
                    auto afmParametersMap = data.getDataMapIgnoringLabels(labels);

                    for(auto& param : afmParametersMap){
                        afmParam_h.push_back(afm->processAFMParameters(param));
                    }

                    //Copy to gpu
                    afmParam_d = afmParam_h;

                    //Check groupslist and afmParameters have the same size
                    if(tip->getNGroups() != afmParam_h.size()){
                        System::log<System::CRITICAL>("[AFMInteractor] (%s) Tip and afmParameters have different sizes (%d != %d)",
                                                       name.c_str(),tip->getNGroups(),afmParam_h.size());
                    }
                    if(sample->getNGroups() != afmParam_h.size()){
                        System::log<System::CRITICAL>("[AFMInteractor] (%s) Sample and afmParameters have different sizes (%d != %d)",
                                                       name.c_str(),sample->getNGroups(),afmParam_h.size());
                    }
                }

                //Set up buffers
                {
                    auto tipParam_h    = tip->getGroupsListInfo(access::location::cpu);
                    auto sampleParam_h = sample->getGroupsListInfo(access::location::cpu);

                    thrust::host_vector<int> id2BufferIndex_h(pd->getNumParticles());
                    thrust::host_vector<int> bufferOffset_h(nAFM+1);

                    //Set id2BufferIndex_h to -1
                    thrust::fill(id2BufferIndex_h.begin(),id2BufferIndex_h.end(),-1);

                    int bufferIndex = 0;
                    for(int i = 0; i < nAFM; ++i){
                        bufferOffset_h[i] = bufferIndex;

                        id2BufferIndex_h[tipParam_h.listStart[i][0]] = bufferIndex;
                        bufferIndex++;
                        for(int j = 0; j < sampleParam_h.n[i]; ++j){
                            id2BufferIndex_h[sampleParam_h.listStart[i][j]] = bufferIndex;
                            bufferIndex++;
                        }
                    }
                    bufferOffset_h[nAFM] = bufferIndex;

                    energyBuffer_d.resize(bufferIndex);
                    forceBuffer_d.resize(bufferIndex);

                    //Copy to gpu
                    id2BufferIndex_d = id2BufferIndex_h;
                    bufferOffset_d   = bufferOffset_h;

                }

                //Set up cub reduction
                {
                    cubReductionTempStorage     = NULL;
                    cubEnergyReductionTempStorageSize = 0;

                    cub::DeviceSegmentedReduce::Reduce(cubReductionTempStorage,
                                                       cubEnergyReductionTempStorageSize,
                                                       energyBuffer_d.begin(),
                                                       getEnergyPermutationIterator(),
                                                       nAFM,
                                                       bufferOffset_d.begin(),
                                                       bufferOffset_d.begin()+1,
                                                       cub::Sum(),
                                                       real(0.0), //Initial value
                                                       0);

                    cubReductionTempStorage     = NULL;
                    cubForceReductionTempStorageSize = 0;

                    cub::DeviceSegmentedReduce::Reduce(cubReductionTempStorage,
                                                       cubForceReductionTempStorageSize,
                                                       forceBuffer_d.begin(),
                                                       getForcePermutationIterator(),
                                                       nAFM,
                                                       bufferOffset_d.begin(),
                                                       bufferOffset_d.begin()+1,
                                                       cub::Sum(),
                                                       make_real4(0.0), //Initial value
                                                       0);

                    cudaMalloc(&cubReductionTempStorage,
                                cubForceReductionTempStorageSize*sizeof(real4));

                }

                cudaDeviceSynchronize();

            }

            ~AFMInteractor(){
                if(cubReductionTempStorage != NULL){
                    cudaFree(cubReductionTempStorage);
                }
            }

            std::shared_ptr<AFMType> getAFM(){
                return afm;
            }
            std::vector<AFMParameters> getAFMParameters(){
                return std::vector<AFMParameters>(afmParam_h.begin(),afmParam_h.end());
            }

            void sum(Computables comp,cudaStream_t st) override {

                auto tipParam    = tip->getGroupsListInfo(access::location::gpu);
                auto sampleParam = sample->getGroupsListInfo(access::location::gpu);

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<AFMType>::value){

                        int Nthreads = THREADS_PER_BLOCK;
                        int Nblocks  = nParticlesSample/Nthreads + ((nParticlesSample%Nthreads)?1:0);
                        if(Nblocks == 0) {Nblocks = 1;} //If there are no particles, we still need to run the kernel
                        // to compute the energy of the tip

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        auto batId = pd->getBatchId(access::location::gpu,access::mode::read).raw();

                        AFMInteractor_ns::transverseAFM
                        <typename AFMType::ComputationalData,
                         typename AFMType::AFMParameters,
                         typename AFMType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(nAFM,
                                                     nParticlesSample,
                                                     tipParam.list,
                                                     sampleParam.list,
                                                     id2index,
                                                     batId,
                                                     afm->getComputationalData(),
                                                     thrust::raw_pointer_cast(afmParam_d.data()),
                                                     thrust::raw_pointer_cast(id2BufferIndex_d.data()),
                                                     thrust::raw_pointer_cast(energyBuffer_d.data()),
                                                     afm->getEnergyTransverser(),
                                                     false);
                        CudaCheckError();

                        cub::DeviceSegmentedReduce::Reduce(cubReductionTempStorage,
                                                           cubEnergyReductionTempStorageSize,
                                                           energyBuffer_d.begin(),
                                                           getEnergyPermutationIterator(),
                                                           nAFM,
                                                           bufferOffset_d.begin(),
                                                           bufferOffset_d.begin()+1,
                                                           cub::Sum(),
                                                           real(0.0), //Initial value
                                                           st);
                        CudaCheckError();
                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[AFMInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }
                }

                if(comp.force == true){

                    if constexpr (has_getForceTransverser<AFMType>::value){

                        int Nthreads = THREADS_PER_BLOCK;
                        int Nblocks  = nParticlesSample/Nthreads + ((nParticlesSample%Nthreads)?1:0);
                        if(Nblocks == 0) {Nblocks = 1;} //If there are no particles, we still need to run the kernel
                        // to compute the force of the tip

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        auto batId = pd->getBatchId(access::location::gpu,access::mode::read).raw();

                        AFMInteractor_ns::transverseAFM
                        <typename AFMType::ComputationalData,
                         typename AFMType::AFMParameters,
                         typename AFMType::ForceTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(nAFM,
                                                     nParticlesSample,
                                                     tipParam.list,
                                                     sampleParam.list,
                                                     id2index,
                                                     batId,
                                                     afm->getComputationalData(),
                                                     thrust::raw_pointer_cast(afmParam_d.data()),
                                                     thrust::raw_pointer_cast(id2BufferIndex_d.data()),
                                                     thrust::raw_pointer_cast(forceBuffer_d.data()),
                                                     afm->getForceTransverser(),
                                                     true);
                        CudaCheckError();

                        cub::DeviceSegmentedReduce::Reduce(cubReductionTempStorage,
                                                           cubForceReductionTempStorageSize,
                                                           forceBuffer_d.begin(),
                                                           getForcePermutationIterator(),
                                                           nAFM,
                                                           bufferOffset_d.begin(),
                                                           bufferOffset_d.begin()+1,
                                                           cub::Sum(),
                                                           make_real4(0.0), //Initial value
                                                           st);
                        CudaCheckError();
                    } else {
                        if(!warningForce){
                            System::log<System::WARNING>("[AFMInteractor] (%s) Requested non-implemented transverser (force)",
                                                         name.c_str());
                            warningForce = true;
                        }
                    }
                }
            }

    };

}}}

#endif
