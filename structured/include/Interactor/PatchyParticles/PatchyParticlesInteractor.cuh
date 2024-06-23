#ifndef __PATCHY_PARTICLES__
#define __PATCHY_PARTICLES__

//TODO Check groups

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

#include "Interactor/Patches/PatchesIncluders.cuh"
#include "Interactor/Patches/GenericPatchesPotentialLoader.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace PatchyParticles_ns{

        __global__ void updatePatchPositions(int numberPatches,
                                             const int*    __restrict__ ppParentId,
                                             const int*    __restrict__ id2index,
                                             const real4*  __restrict__ pos,
                                             const real4*  __restrict__ dir,
                                             const real4*  __restrict__ ppPatchPos,
                                                   real4*  __restrict__ ppPatchVector,
                                                   real4*  __restrict__ ppPos,
                                                   int*    __restrict__ ppParentIndex){

            int patchIndex = blockIdx.x*blockDim.x + threadIdx.x;

            if(patchIndex>=numberPatches){return;}

            int parentIndex = id2index[ppParentId[patchIndex]];

            Quat q = dir[parentIndex];

            const real4 localCoord = ppPatchPos[patchIndex];
            const real3 localPos = rotateVector(q,make_real3(localCoord));

            ppPatchVector[patchIndex].x=localPos.x;
            ppPatchVector[patchIndex].y=localPos.y;
            ppPatchVector[patchIndex].z=localPos.z;

            ppPos[patchIndex].x=pos[parentIndex].x+localPos.x;
            ppPos[patchIndex].y=pos[parentIndex].y+localPos.y;
            ppPos[patchIndex].z=pos[parentIndex].z+localPos.z;

            ppParentIndex[patchIndex] = parentIndex;

        }

        __global__ void accumulateRealAndZero(int numberParticles,
                                              ParticleGroup::IndexIterator groupIndex,
                                              const int*  __restrict__ partId,
                                              const int*  __restrict__ id2index_patches,
                                                    int** __restrict__ listStart,
                                                    int*  __restrict__ n,
                                                    real*  __restrict__ quantityParticle,
                                                    real*  __restrict__ quantityPatches){

            int localIndex = blockIdx.x*blockDim.x + threadIdx.x;

            if(localIndex>=numberParticles){return;}

            const int partIndex = groupIndex[localIndex];
            const int id        = partId[partIndex];

            const int* start = listStart[id];
            real acc = real(0.0);
            for(int i=0;i<n[id];i++){
                const int patchIndex = id2index_patches[start[i]];
                acc+=quantityPatches[patchIndex];
                quantityPatches[patchIndex]=real(0.0);
            }
            quantityParticle[partIndex]+=acc;
        }

        __global__ void accumulateReal4Real4AndZero(int numberParticles,
                                                    ParticleGroup::IndexIterator groupIndex,
                                                    const int*  __restrict__ partId,
                                                    const int*  __restrict__ id2index_patches,
                                                          int** __restrict__ listStart,
                                                          int*  __restrict__ n,
                                                          real4*  __restrict__ quantityParticle1,
                                                          real4*  __restrict__ quantityPatches1,
                                                          real4*  __restrict__ quantityParticle2,
                                                          real4*  __restrict__ quantityPatches2){

            int localIndex = blockIdx.x*blockDim.x + threadIdx.x;

            if(localIndex>=numberParticles){return;}

            const int partIndex = groupIndex[localIndex];
            const int id        = partId[partIndex];

            const int* start = listStart[id];
            real4 acc1 = make_real4(0.0);
            real4 acc2 = make_real4(0.0);
            for(int i=0;i<n[id];i++){
                const int patchIndex = id2index_patches[start[i]];
                acc1+=quantityPatches1[patchIndex];
                acc2+=quantityPatches2[patchIndex];
                quantityPatches1[patchIndex]=make_real4(0.0);
                quantityPatches2[patchIndex]=make_real4(0.0);
            }
            quantityParticle1[partIndex]+=acc1;
            quantityParticle2[partIndex]+=acc2;
            //printf("%i f: %f %f %f %f    t: %f %f %f %f\n",id,acc1.x,acc1.y,acc1.z,acc1.w,acc2.x,acc2.y,acc2.z,acc2.w);
        }
    }

    template<int THREADS_PER_BLOCK = 256>
    class PatchyParticles_: public Interactor{

        protected:

            //Particles

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;

            //Patches

            std::shared_ptr<GlobalData>           patchesGd;
            std::shared_ptr<ExtendedParticleData> patchesPd;

            //////////////////////////////

            std::shared_ptr<InputEntryManager>    patchForceFieldInfo;

            //////////////////////////////

            std::map<std::string, std::shared_ptr<ParticleGroup>> patchesGroups;
            std::map<std::string, std::shared_ptr<VerletConditionalListSetBase>> patchesVConListSet;
            std::map<std::string, std::shared_ptr<Interactor>> patchInteractors;

            std::shared_ptr<groupsList>  particleId2patchesId;

            void accumulate(uammd::Interactor::Computables comp,cudaStream_t st) {

                auto patchesList = particleId2patchesId->getGroupsListInfo(access::location::gpu);

                int numberParticles     = this->pg->getNumberParticles();
                auto groupIndexIterator = this->pg->getIndexIterator(access::location::gpu);

                ////////////////////////////////////////////////////////////////

                auto id = pd->getId(access::location::gpu, access::mode::read);
                auto id2index_patches = patchesPd->getIdOrderedIndices(access::location::gpu);

                int** listStart = patchesList.listStart;
                int*  n = patchesList.n;

                /////////////////////////////////////////////////////////////////////////////

                int Nthreads=THREADS_PER_BLOCK;
                int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                /////////////////////////////////////////////////////////////////////////////


                if(comp.energy == true){

                    auto energyPart  = pd->getEnergy(access::location::gpu, access::mode::readwrite);
                    auto energyPatch = patchesPd->getEnergy(access::location::gpu, access::mode::read);

                    PatchyParticles_ns::accumulateRealAndZero
                    <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                 groupIndexIterator,
                                                 id.raw(),
                                                 id2index_patches,
                                                 listStart,
                                                 n,
                                                 energyPart.raw(),
                                                 energyPatch.raw());
                    CudaCheckError();
                }

                if(comp.force == true){

                    auto forcePart  = pd->getForce(access::location::gpu, access::mode::readwrite);
                    auto forcePatch = patchesPd->getForce(access::location::gpu, access::mode::read);

                    auto torquePart  = pd->getTorque(access::location::gpu, access::mode::readwrite);
                    auto torquePatch = patchesPd->getTorque(access::location::gpu, access::mode::read);

                    PatchyParticles_ns::accumulateReal4Real4AndZero
                    <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                 groupIndexIterator,
                                                 id.raw(),
                                                 id2index_patches,
                                                 listStart,
                                                 n,
                                                 forcePart.raw(),
                                                 forcePatch.raw(),
                                                 torquePart.raw(),
                                                 torquePatch.raw());
                    CudaCheckError();
                }
            }


        public:

            PatchyParticles_(std::shared_ptr<GlobalData>    gd,
                             std::shared_ptr<ParticleGroup> pg,
                             std::vector<std::string>     path,
                             std::string name):Interactor(pg,"PatchyParticlesInteractor: \"" +name+"\""),
                                               gd(gd){

                pd  = getExtendedParticleData(this->pg);

                //Check if pg is all (temporal)
                if(pg->getNumberParticles() != pd->getNumParticles()){
                    System::log<System::CRITICAL>("[PatchyParticles] (%s) The version of patchy particles for a group"
                                                  " different than particle data has not been implemented yet.",name.c_str());
                }

                //Load patches
                //Load patches global data
                std::vector<std::string> pathPatchTypes = path;
                pathPatchTypes.push_back("patchesGlobal");
                this->patchesGd = std::make_shared<GlobalData>(gd->getGlobalDataBase(),pathPatchTypes);

                //Load patches positions
                std::vector<std::string> pathPatchPositions = path;
                pathPatchPositions.push_back("patchesState");
                this->patchesPd = std::make_shared<ExtendedParticleData>(getExtendedSystem(this->sys),pathPatchPositions);

                ///////////////////////////////////
                //     Load patches topology     //
                ///////////////////////////////////

                std::vector<std::string> pathPatchTopology = path;
                pathPatchTopology.push_back("patchesTopology");

                //Load structure
                std::vector<std::string> availStructureLabels = {"id","type","parentId","batchId"};

                std::vector<std::string> pathPatchStructure = pathPatchTopology;
                pathPatchStructure.push_back("structure");

                DataEntry structureData = (getExtendedSystem(this->sys))->
                                           getInput()->
                                           getDataEntry(pathPatchStructure);

                //Check all labels of structure data are available
                for(auto label : structureData.getLabels()){
                    if(std::find(availStructureLabels.begin(),availStructureLabels.end(),label) == availStructureLabels.end()){
                        System::log<System::CRITICAL>("[PatchyParticles] (%s) The label \"%s\" is not available in the structure data.",
                                                      name.c_str(),label.c_str());
                    }
                }

                std::vector<int>         id       = structureData.getData<int>("id");
                std::vector<std::string> type     = structureData.getData<std::string>("type");
                std::vector<int>         parentId = structureData.getData<int>("parentId");

                bool isBatchIdGiven = structureData.isDataAdded("batchId");
                std::vector<int> batchId;
                if(isBatchIdGiven){
                    batchId = structureData.getData<int>("batchId");
                }

                auto patchTypeParamHandler = this->patchesGd->getTypes();

                int N = id.size();

                //Check N == patchesPd->getNumParticles()
                if(N != patchesPd->getNumParticles()){
                    System::log<System::CRITICAL>("[PatchyParticles] (%s) The number of patches in the structure entry"
                                                  " is different than the number of patches in the patches coordinates.",name.c_str());
                }

                {
                    auto ppId       = patchesPd->getId(access::location::cpu, access::mode::read);
                    auto ppPos      = patchesPd->getPos(access::location::cpu, access::mode::write);
                    auto ppParentId = patchesPd->getModelId(access::location::cpu, access::mode::write);
                    auto ppBatId    = patchesPd->getBatchId(access::location::cpu, access::mode::write);

                    //The batch of each patch is equal to the batch id of the parent particle
                    auto id2index   = pd->getIdOrderedIndices(access::location::cpu);
                    auto batId      = pd->getBatchId(access::location::cpu, access::mode::read);

                    for(int index=0;index<N;index++){
                        if(ppId[index] != id[index]){
                            System::log<System::CRITICAL>("[PatchyParticles] (%s) The patch id %d does not match with the id %d in the structure entry.",
                                                          name.c_str(),ppId[index],id[index]);
                        } else {
                            ppPos[index].w    = int(patchTypeParamHandler->getTypeId(type[index]));
                            ppParentId[index] = parentId[index];
                            ppBatId[index]    = batId[id2index[parentId[index]]];

                            if(isBatchIdGiven){
                                if(ppBatId[index] != batchId[index]){
                                    System::log<System::CRITICAL>("[PatchyParticles] (%s) The batch id %d of the patch %d "
                                                                  "( equal to the batch id of the parent particle %d )"
                                                                  " does not match with the batch id %d in the structure entry.",
                                                                  name.c_str(),ppBatId[index],ppId[index],ppParentId[index],batchId[index]);
                                }
                            }
                        }
                    }
                }

                {
                    System::log<System::MESSAGE>("[PatchyParticles] (%s) Loading patches types.",name.c_str());
                    patchTypeParamHandler->loadTypesIntoParticleData(patchesPd);
                }

                //Load force field
                std::vector<std::string> pathPatchForceField = pathPatchTopology;
                pathPatchForceField.push_back("forceField");
                patchForceFieldInfo =
                std::make_shared<InputEntryManager>(getExtendedSystem(this->sys),
                                                    pathPatchForceField);

                //Load groups
                patchesGroups = GroupUtils::loadGroupsListFromInputEntries(getExtendedSystem(this->sys),
                                                                           patchesGd,patchesPd,
                                                                           patchForceFieldInfo);

                //Load neig lists
                patchesVConListSet = VerletConditionalListSetUtils::loadVerletConditionalListSetsFromInputEntries(getExtendedSystem(this->sys),
                                                                                                                  patchesGd,patchesGroups,
                                                                                                                  patchForceFieldInfo);
                //Load patches interactors

                for(auto& entry : patchForceFieldInfo->getEntriesInfo()){

                    if(Potentials::GenericPatchyParticlesLoader::isPatchyParticlesInteractorAvailable(getExtendedSystem(this->sys),entry.second.path)){
                      std::shared_ptr<typename uammd::Interactor> inter = Potentials::GenericPatchyParticlesLoader::loadGenericPatchyParticles(getExtendedSystem(this->sys),
                                                                                                                                               gd,pg,
                                                                                                                                               patchesGd,patchesGroups,
                                                                                                                                               patchesVConListSet,entry.second.path);

                      if(patchInteractors.count(entry.second.name) == 0){
                            patchInteractors[entry.second.name] = inter;
                            entry.second.used = true;
                        } else {
                            System::log<System::CRITICAL>("[PatchyParticles] (%s) Error loading interactors,"
                                                          "interactor \"%s\" has already been added.",name.c_str(),entry.second.name.c_str());
                        }
                    }
                }

                //Print information about interactors
                for(auto i : patchInteractors){
                    System::log<System::MESSAGE>("[PatchyParticles] (%s) Added interactor \"%s\".",
                                                 name.c_str(),i.first.c_str());
                }

                //Check used entries

                patchForceFieldInfo->checkEntriesUsed();

                ///////////////////////////////////

                //Copy positions to patchePositions
                {
                    auto ppPos    = patchesPd->getPos(access::location::cpu, access::mode::read);
                    auto ppPatPos = patchesPd->getPatchPos(access::location::cpu, access::mode::write);

                    thrust::copy(ppPos.begin(), ppPos.end(), ppPatPos.begin());
                }

                //Create mol to patches map
                std::vector<int>              ids(pd->getNumParticles());
                std::vector<std::vector<int>> patches(pd->getNumParticles());

                {
                    auto pId       = pd->getId(access::location::cpu, access::mode::read);

                    auto ppId      = patchesPd->getId(access::location::cpu, access::mode::read);
                    auto ppParentId = patchesPd->getModelId(access::location::cpu, access::mode::read);

                    for(int i=0;i<pd->getNumParticles();i++){
                        ids[i]=pId[i];
                        if(i != pId[i]){
                            System::log<System::CRITICAL>("[PatchyParticles] (%s) Index and particle id have to match",name.c_str());
                        }
                    }

                    for(int i=0;i<patchesPd->getNumParticles();i++){
                        int particleId = ppParentId[i];
                        patches[particleId].push_back(ppId[i]);
                    }
                }

                particleId2patchesId = std::make_shared<groupsList>(ids,patches);

                //for(int i=0;i<pd->getNumParticles();i++){
                //    std::cout << ids[i] << ": ";
                //    for(int pid : patches[i]){
                //        std::cout << pid << " ";
                //    }
                //    std::cout << std::endl;
                //}

            }

            void updatePatchyParticles(cudaStream_t st){

                auto ppParentId = patchesPd->getModelId(access::location::gpu, access::mode::read);
                auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                auto pos      = pd->getPos(access::location::gpu, access::mode::read);
                auto dir      = pd->getDir(access::location::gpu, access::mode::read);

                auto ppPatchPos     = patchesPd->getPatchPos(access::location::gpu, access::mode::read);
                auto ppPatchVector  = patchesPd->getPatchVector(access::location::gpu, access::mode::write);
                auto ppPos          = patchesPd->getPos(access::location::gpu, access::mode::write);

                auto ppParentIndex  = patchesPd->getParentIndex(access::location::gpu, access::mode::write);

                /////////////////////////////////////////////////////////////////////////////

                int numberPatches = patchesPd->getNumParticles();

                int Nthreads=THREADS_PER_BLOCK;
                int Nblocks=numberPatches/Nthreads + ((numberPatches%Nthreads)?1:0);

                PatchyParticles_ns::updatePatchPositions
                <<<Nblocks,Nthreads,0, st>>>(numberPatches,
                                             ppParentId.raw(),
                                             id2index,
                                             pos.raw(),
                                             dir.raw(),
                                             ppPatchPos.raw(),
                                             ppPatchVector.raw(),
                                             ppPos.raw(),
                                             ppParentIndex.raw());
                CudaCheckError();
            }


            std::shared_ptr<GlobalData>                         getPatchesGlobalData()  {return patchesGd;}
            std::shared_ptr<ExtendedParticleData>               getPatchesParticleData(){return patchesPd;}
            std::map<std::string, std::shared_ptr<Interactor>>  getPatchesInteractors() {return patchInteractors;}

            InputEntryManager::entryInfo                        getPatchesInteractorInfo(std::string name){
                if(patchInteractors.count(name) == 0){
                    System::log<System::CRITICAL>("[PatchyParticles] Error getting interactor data \"%s\", interactor has not been added.",name.c_str());
                }
                return patchForceFieldInfo->getEntryInfo(name);
            }



            virtual void sum(uammd::Interactor::Computables comp,cudaStream_t st) override {
                this->updatePatchyParticles(st);
                for(auto &interactor: this->patchInteractors){
                    interactor.second->sum(comp,st);
                }
                this->accumulate(comp,st);
            }
    };

    namespace DynamicallyBondedPatchyParticles_ns{

        __global__ void updatePatchState(int numberPatches,
                                         const int*   __restrict__ id2index,
                                         const int*   __restrict__ patchesId,
                                         int4*        __restrict__ ppState,
                                         real4*       __restrict__ ppStateInfo,
                                         real energyThreshold){

            int index = blockIdx.x*blockDim.x + threadIdx.x;

            if(index>=numberPatches){return;}

            const real energy = ppStateInfo[index].x;

            if(energy<energyThreshold){
                // Tentative bond
                // if energy is less than threshold,
                // then it is ensured that ppState[index].y != -1
                const int tentativeBondId    = ppState[index].y;
                const int tentativeBondIndex = id2index[tentativeBondId];

                const int tentBondTentBond = ppState[tentativeBondIndex].y;

                // Check if the tentativeBond of tentativeBond is the current patch id
                if(tentBondTentBond == patchesId[index]){
                    ppState[index].x = tentativeBondId;
                } else {
                    ppState[index].x = int(-1);
                }
            } else {
                ppState[index].x = int(-1);
            }

            ppStateInfo[index].x = energyThreshold;

        }

        __global__ void resetPatchState(int numberPatches,
                                        int4* __restrict__ ppState){

            int index = blockIdx.x*blockDim.x + threadIdx.x;

            if(index>=numberPatches){return;}

            ppState[index].y = int(-1);

        }


    }

    template<int THREADS_PER_BLOCK = 256>
    class DynamicallyBondedPatchyParticles_: public PatchyParticles_<THREADS_PER_BLOCK>{

        private:

            using Base = PatchyParticles_<THREADS_PER_BLOCK>;

            real energyThreshold;

        public:

            DynamicallyBondedPatchyParticles_(std::shared_ptr<GlobalData>    gd,
                                              std::shared_ptr<ParticleGroup> pg,
                                              std::vector<std::string>     path,
                                              std::string name):Base(gd,pg,path,name){

                energyThreshold = this->patchesGd->getFundamental()->getEnergyThreshold();

                System::log<System::MESSAGE>("[DynamicallyBondedPatchyParticles] (%s) energyThreshold = %f",
                                             name.c_str(),energyThreshold);

                /////////////////////////////////////////////////////////////////////////////////////

                auto ppState     = this->patchesPd->getState(access::location::cpu, access::mode::write);
                auto ppStateInfo = this->patchesPd->getStateInfo(access::location::cpu, access::mode::write);

                int4 initInt4 = {-1,-1,0,0};
                std::fill(ppState.begin(),ppState.end(),initInt4);
                std::fill(ppStateInfo.begin(),ppStateInfo.end(),make_real4(0.0));
            }

            void updatePatchyParticles(cudaStream_t st){
                Base::updatePatchyParticles(st);

                {

                    int numberPatches = this->patchesPd->getNumParticles();

                    auto id2index    = this->patchesPd->getIdOrderedIndices(access::location::gpu);
                    auto ppPatchesId = this->patchesPd->getId(access::location::gpu, access::mode::read);
                    auto ppState     = this->patchesPd->getState(access::location::gpu, access::mode::readwrite);
                    auto ppStateInfo = this->patchesPd->getStateInfo(access::location::gpu, access::mode::readwrite);

                    int Nthreads=THREADS_PER_BLOCK;
                    int Nblocks=numberPatches/Nthreads + ((numberPatches%Nthreads)?1:0);

                    DynamicallyBondedPatchyParticles_ns::updatePatchState
                    <<<Nblocks,Nthreads,0, st>>>(numberPatches,
                                                 id2index,
                                                 ppPatchesId.raw(),
                                                 ppState.raw(),
                                                 ppStateInfo.raw(),
                                                 energyThreshold);
                    CudaCheckError();
                }

                {
                    int numberPatches = this->patchesPd->getNumParticles();

                    auto ppState = this->patchesPd->getState(access::location::gpu, access::mode::readwrite);

                    int Nthreads=THREADS_PER_BLOCK;
                    int Nblocks=numberPatches/Nthreads + ((numberPatches%Nthreads)?1:0);

                    DynamicallyBondedPatchyParticles_ns::resetPatchState
                    <<<Nblocks,Nthreads,0, st>>>(numberPatches,
                                                 ppState.raw());
                    CudaCheckError();
                }
            }

            void sum(uammd::Interactor::Computables comp,cudaStream_t st) override {
                this->updatePatchyParticles(st);
                for(auto &interactor: this->patchInteractors){
                    interactor.second->sum(comp,st);
                }
                this->accumulate(comp,st);
            }

    };

}}}

#endif
