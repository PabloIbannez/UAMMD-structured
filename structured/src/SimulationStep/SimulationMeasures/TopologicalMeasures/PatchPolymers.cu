#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

#include <thrust/iterator/transform_output_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Definitions/Matrices.cuh"

#include "Interactor/PatchyParticles/PatchyParticles/PatchyParticles.cu"
#include "Interactor/PatchyParticles/PatchyParticles/DynamicallyBondedPatchyParticles.cu"

#include "Interactor/Patches/SurfacePatches/Linker.cu"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

        namespace PatchPolymers_ns{

                enum eventType {Pp=0, Pn=1, Dp=2, Dn=3, S=4, B=5};
                //Pp: polimerization positive direction
                //Pn: polimerization negative direction
                //Dp: depolimerization positive direction
                //Dn: depolimerization negative direction
                //S:  from bulk to surface
                //B:  from surface to bulk

                //Negative direction events (Pn and Dn) are ignored.
                //The reason is that the polymerization and depolimerization in the negative direction
                //can be inferred from the positive direction events.


                #pragma pack(push, 1) //Remove padding bytes
                struct event{
                        ullint    step; //step in which the event occurs
                        int       id;   //id of the particle (monomer)
                        eventType type; //type of event
                        int       info; //info about the event
                };
                #pragma pack(pop)

                __device__ inline void recordPolimerizationEvents(int id,
                                ullint step,
                                int Nmonomers,
                                ullint bufferSize,
                                const int* id2indexMonomers,
                                const int* id2indexPatches,
                                const int2 * monomer2patches,
                                const int  * patch2monomer,
                                int4 * previousState,
                                const int4 * currentState,
                                event* eventBuffer,
                                int*   nEventsPerMonomer){

                        const int endPatchId = monomer2patches[id].x;

                        //Check endPatch
                        //const int prevEndState = previousState[endPatchId].x;
                        const int currEndState = currentState[id2indexPatches[endPatchId]].x;

                        //if(prevEndState != currEndState){
                        //	//Event detected
                        //	if(prevEndState == int(-1)){
                        //		//Polimerization
                        //		event newEvent;

                        //		newEvent.step = step;
                        //		newEvent.id   = id;
                        //		newEvent.type = Pn;
                        //		newEvent.info = patch2monomer[currEndState];

                        //		UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                        //													 bufferSize, Nmonomers,
                        //													 nEventsPerMonomer[id], id,
                        //													 newEvent);
                        //		nEventsPerMonomer[id]++;
                        //	}

                        //	if(currEndState == int(-1)){
                        //		//Depolimerization
                        //		event newEvent;

                        //		newEvent.step = step;
                        //		newEvent.id   = id;
                        //		newEvent.type = Dn;
                        //		newEvent.info = patch2monomer[prevEndState];

                        //		UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                        //													 bufferSize, Nmonomers,
                        //													 nEventsPerMonomer[id], id,
                        //													 newEvent);
                        //		nEventsPerMonomer[id]++;
                        //	}
                        //}

                        const int begPatchId = monomer2patches[id].y;

                        //Check begPatch
                        const int prevBegState = previousState[begPatchId].x;
                        const int currBegState = currentState[id2indexPatches[begPatchId]].x;

                        if(prevBegState != currBegState){
                                //Event detected
                                if(prevBegState == int(-1)){
                                        //Polimerization
                                        event newEvent;

                                        newEvent.step = step;
                                        newEvent.id   = id;
                                        newEvent.type = Pp;
                                        newEvent.info = patch2monomer[currBegState];

                                        UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                                                        bufferSize, Nmonomers,
                                                        nEventsPerMonomer[id], id,
                                                        newEvent);
                                        nEventsPerMonomer[id]++;
                                }

                                if(currBegState == int(-1)){
                                        //Depolimerization
                                        event newEvent;

                                        newEvent.step = step;
                                        newEvent.id   = id;
                                        newEvent.type = Dp;
                                        newEvent.info = patch2monomer[prevBegState];

                                        UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                                                        bufferSize, Nmonomers,
                                                        nEventsPerMonomer[id], id,
                                                        newEvent);
                                        nEventsPerMonomer[id]++;
                                }
                        }

                        //Update

                        previousState[endPatchId].x = currEndState;
                        previousState[begPatchId].x = currBegState;
                }

                __device__ inline void recordSurfaceEvents(int id,
                                ullint step,
                                int Nmonomers,
                                ullint bufferSize,
                                const int* id2indexMonomers,
                                const int* id2indexSurfPatches,
                                const int  * monomer2surfPatches,
                                real * previousSurfaceEnergy,
                                typename Potentials::SurfacePatches::Linker::EnergyTransverser transverser,
                                typename Potentials::SurfacePatches::Linker::ComputationalData computationalData,
                                real surfaceEnergyThreshold,
                                event* eventBuffer,
                                int*   nEventsPerMonomer){

                        const int surfacePatchId    = monomer2surfPatches[id];
                        const int indexSurfacePatch = id2indexSurfPatches[surfacePatchId];

                        //Check surface events
                        const real prevSurfaceEnergy = previousSurfaceEnergy[id];
                        const real currSurfaceEnergy = transverser.compute(indexSurfacePatch,computationalData);

                        if(prevSurfaceEnergy >= surfaceEnergyThreshold and
                                        currSurfaceEnergy <  surfaceEnergyThreshold){
                                //Bulk to surface
                                event newEvent;

                                newEvent.step = step;
                                newEvent.id   = id;
                                newEvent.type = S;
                                newEvent.info = 0;

                                UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                                                bufferSize, Nmonomers,
                                                nEventsPerMonomer[id], id,
                                                newEvent);
                                nEventsPerMonomer[id]++;
                        }

                        if(prevSurfaceEnergy < surfaceEnergyThreshold and
                                        currSurfaceEnergy >= surfaceEnergyThreshold){
                                //Surface to bulk
                                event newEvent;

                                newEvent.step = step;
                                newEvent.id   = id;
                                newEvent.type = B;
                                newEvent.info = 0;

                                UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                                                bufferSize, Nmonomers,
                                                nEventsPerMonomer[id], id,
                                                newEvent);
                                nEventsPerMonomer[id]++;
                        }

                        //Update

                        previousSurfaceEnergy[id] = currSurfaceEnergy;

                }

                __global__ void recordEvents(const ullint step,
                                const int Nmonomers,
                                const ullint bufferSize,
                                const int* id2indexMonomers,
                                const int* id2indexPatches,
                                const int2 * __restrict__ monomer2patches,
                                const int  * __restrict__ patch2monomer,
                                int4 * __restrict__ previousState,
                                const int4 * __restrict__ currentState,
                                event* __restrict__ eventBuffer,
                                int*   __restrict__ nEventsPerMonomer,
                                int*                maxEvents){

                        //One thread per monomer
                        const int id = blockIdx.x*blockDim.x + threadIdx.x;

                        if(id >= Nmonomers) return;

                        recordPolimerizationEvents(id,
                                        step,
                                        Nmonomers,
                                        bufferSize,
                                        id2indexMonomers,id2indexPatches,
                                        monomer2patches,patch2monomer,
                                        previousState,
                                        currentState,
                                        eventBuffer,nEventsPerMonomer);

                        atomicMax(maxEvents, nEventsPerMonomer[id]);
                }

                __global__ void recordEventsSurf(const ullint step,
                                const int Nmonomers,
                                const ullint bufferSize,
                                const int* id2indexMonomers,
                                const int* id2indexPatches,
                                const int* id2indexSurfPatches,
                                const int2 * __restrict__ monomer2patches,
                                const int  * __restrict__ patch2monomer,
                                const int  * __restrict__ monomer2surfPatches,
                                int4 * __restrict__ previousState,
                                const int4 * __restrict__ currentState,
                                real*  __restrict__ previousSurfaceEnergy,
                                typename Potentials::SurfacePatches::Linker::EnergyTransverser transverser,
                                typename Potentials::SurfacePatches::Linker::ComputationalData computationalData,
                                const real surfaceEnergyThreshold,
                                event* __restrict__ eventBuffer,
                                int*   __restrict__ nEventsPerMonomer,
                                int*                maxEvents){

                        //One thread per monomer
                        const int id = blockIdx.x*blockDim.x + threadIdx.x;

                        if(id >= Nmonomers) return;

                        recordPolimerizationEvents(id,
                                        step,
                                        Nmonomers,
                                        bufferSize,
                                        id2indexMonomers,id2indexPatches,
                                        monomer2patches,patch2monomer,
                                        previousState,
                                        currentState,
                                        eventBuffer,nEventsPerMonomer);

                        recordSurfaceEvents(id,
                                        step,
                                        Nmonomers,
                                        bufferSize,
                                        id2indexMonomers,id2indexSurfPatches,
                                        monomer2surfPatches,
                                        previousSurfaceEnergy,
                                        transverser,
                                        computationalData,
                                        surfaceEnergyThreshold,
                                        eventBuffer,
                                        nEventsPerMonomer);

                        atomicMax(maxEvents, nEventsPerMonomer[id]);
                }
        }

        class PatchPolymers : public SimulationStepBase{

                private:

                        int THREADS_PER_BLOCK = 512;

                        std::string   outputFilePath;
                        std::ofstream outputFile;

                        bool isPolymer = false;
                        std::shared_ptr<Interactor::PatchyParticles::DynamicallyBondedPatchyParticles> polymerPatchyParticles;

                        bool isSurface = false;
                        std::shared_ptr<Interactor::PatchyParticles::PatchyParticles> surfacePatchyParticles;
                        std::shared_ptr<Potentials::SurfacePatches::Linker>           linker;

                        std::shared_ptr<ParticleGroup> polymerGroup; //Group of patches
                        std::shared_ptr<ParticleGroup> surfaceGroup; //Group of patches

                        std::vector<PatchPolymers_ns::event> eventWriteBuffer;

                        //Buffering

                        static bool areParametersSet;
                        static bool isBufferSet;

                        static bool firstStep;

                        static std::string startTypeName;
                        static std::string endTypeName;
                        static std::string linkerTypeName;

                        static int startType;
                        static int endType;
                        static int linkerType;

                        static real surfaceEnergyThreshold;

                        static thrust::device_vector<int2> monomer2patches;
                        static thrust::device_vector<int>  patch2monomer;

                        static thrust::device_vector<int>  monomer2surfPatches;

                        static ullint bufferSize;

                        static ullint lastBufferStep;
                        static ullint lastUpdateStep;

                        static thrust::host_vector<PatchPolymers_ns::event>   eventBufferCPU;
                        static thrust::device_vector<PatchPolymers_ns::event> eventBufferGPU;

                        static thrust::host_vector<int>   nEventsPerMonomerCPU;
                        static thrust::device_vector<int> nEventsPerMonomerGPU;

                        static int*                        maxEventsCPU;
                        static thrust::device_vector<int>  maxEventsGPU;

                        static thrust::device_vector<int4> previousState;
                        static thrust::device_vector<real> previousSurfaceEnergy;

                        static thrust::host_vector<int>                       nEventsPerMonomerRead;
                        static thrust::host_vector<PatchPolymers_ns::event>   eventBufferRead;

                        void writeEvents(){

                                //Get the ids of the particles in the group
                                auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                access::location::cpu);

                                eventWriteBuffer.resize(0);

                                for(int i=0; i <this->pg->getNumberParticles(); i++){
                                        int id = pgIds[i];
                                        for(int n=0; n < nEventsPerMonomerRead[id]; n++){
                                                PatchPolymers_ns::event e = UAMMD_GET_2D_ROW_MAJOR(eventBufferRead,
                                                                bufferSize+6, this->pd->getNumParticles(),
                                                                n, id);
                                                eventWriteBuffer.push_back(e);
                                        }
                                }

                                //Sort the events by step and id
                                std::sort(eventWriteBuffer.begin(), eventWriteBuffer.end(),
                                                [](const PatchPolymers_ns::event& a, const PatchPolymers_ns::event& b){
                                                if(a.step == b.step) return a.id < b.id;
                                                else return a.step < b.step;
                                                });

                                //Write the events to the file in binary format
                                for(auto e : eventWriteBuffer){
                                        outputFile.write((char*)&e, sizeof(PatchPolymers_ns::event));
                                }
                        }

                        void resetBuffer(cudaStream_t st){

                                eventBufferRead       = eventBufferGPU;
                                nEventsPerMonomerRead = nEventsPerMonomerGPU;

                                /////////////////////////////////////////////

                                maxEventsCPU[0] = 0;
                                cudaMemcpyAsync(thrust::raw_pointer_cast(maxEventsGPU.data()), maxEventsCPU,
                                                sizeof(int), cudaMemcpyHostToDevice,st);

                                std::fill(nEventsPerMonomerCPU.begin(),nEventsPerMonomerCPU.end(),0);
                                nEventsPerMonomerGPU = nEventsPerMonomerCPU;

                                cudaStreamSynchronize(st);

                        }

                public:

                        PatchPolymers(std::shared_ptr<ParticleGroup>             pg,
                                        std::shared_ptr<IntegratorManager>  integrator,
                                        std::shared_ptr<ForceField>                 ff,
                                        DataEntry& data,
                                        std::string name):SimulationStepBase(pg,integrator,ff,data,name){

                                //Read the output file path
                                outputFilePath = data.getParameter<std::string>("outputFilePath");

                                if(!areParametersSet){
                                        startTypeName = data.getParameter<std::string>("startType","S");
                                        endTypeName   = data.getParameter<std::string>("endType"  ,"E");
                                        linkerTypeName= data.getParameter<std::string>("linkerType","L");

                                        surfaceEnergyThreshold = data.getParameter<real>("surfaceEnergyThreshold",0.0);

                                        areParametersSet = true;
                                } else {
                                        if(startTypeName != data.getParameter<std::string>("startType","S")){
                                                System::log<System::CRITICAL>("[PatchPolymers] startType is not the same for all PatchPolymers instances");
                                        }
                                        if(endTypeName != data.getParameter<std::string>("endType","E")){
                                                System::log<System::CRITICAL>("[PatchPolymers] endType is not the same for all PatchPolymers instances");
                                        }
                                        if(linkerTypeName != data.getParameter<std::string>("linkerType","L")){
                                                System::log<System::CRITICAL>("[PatchPolymers] linkerType is not the same for all PatchPolymers instances");
                                        }
                                        if(surfaceEnergyThreshold != data.getParameter<real>("surfaceEnergyThreshold")){
                                                System::log<System::CRITICAL>("[PatchPolymers] surfaceEnergyThreshold is not the same for all PatchPolymers instances");
                                        }
                                }

                                ullint newBufferSize = data.getParameter<ullint>("bufferSize");
                                if(newBufferSize != bufferSize and bufferSize != 0){
                                        //Error
                                        System::log<System::CRITICAL>("[PatchPolymers] Buffer step is not the same for all the instances of the class");
                                } else {
                                        bufferSize = newBufferSize;
                                }

                        }

                        //We have to override the init function.
                        void init(cudaStream_t st) override{

                                bool isFileEmpty = Backup::openFile(this->sys, outputFilePath, outputFile, true);

                                if(not isFileEmpty){
                                        System::log<System::MESSAGE>("[PatchPolymers] Output file (%s) is not empty. Appending data.", outputFilePath.c_str());
                                }

                                //Look for the patchy particles
                                std::map<std::string,std::shared_ptr<typename uammd::Interactor>> patchyParticles = this->topology->getInteractorsByClass("PatchyParticles");

                                //Check if patchyParticles is empty
                                if(patchyParticles.empty()){
                                        System::log<System::CRITICAL>("[PatchPolymers] There are no patchy particles in the system.");
                                }

                                //Iterate over the patchy particles, and detect if patchy and surface particles are present
                                for(auto& pp : patchyParticles){
                                        auto patchyParticlesInfo = this->topology->getInteractorInfo(pp.first);
                                        if(patchyParticlesInfo.entrySubType == "DynamicallyBondedPatchyParticles"){
                                                if(isPolymer){
                                                        System::log<System::CRITICAL>("[PatchPolymers] There are more than one DynamicallyBondedPatchyParticles in the system.");
                                                } else {
                                                        isPolymer = true;
                                                        polymerPatchyParticles = std::dynamic_pointer_cast<Interactor::PatchyParticles::DynamicallyBondedPatchyParticles>(pp.second);
                                                        System::log<System::MESSAGE>("[PatchPolymers] Found DynamicallyBondedPatchyParticles (%s).", pp.first.c_str());
                                                }
                                        } else if(patchyParticlesInfo.entrySubType == "PatchyParticles"){
                                                if(isSurface){
                                                        System::log<System::CRITICAL>("[PatchPolymers] There are more than one PatchyParticles in the system.");
                                                } else {
                                                        isSurface = true;
                                                        surfacePatchyParticles = std::dynamic_pointer_cast<Interactor::PatchyParticles::PatchyParticles>(pp.second);
                                                        System::log<System::MESSAGE>("[PatchPolymers] Found PatchyParticles (%s).", pp.first.c_str());
                                                }
                                        }
                                }

                                if(!isPolymer){
                                        System::log<System::CRITICAL>("[PatchPolymers] There are no DynamicallyBondedPatchyParticles in the system.");
                                } else {

                                        //Set type ids
                                        auto types = polymerPatchyParticles->getPatchesGlobalData()->getTypes();

                                        startType = types->getTypeId(startTypeName);
                                        endType   = types->getTypeId(endTypeName);

                    //Check polymer has the correct type
                    for(auto inte : polymerPatchyParticles->getPatchesInteractors()){
                            if( polymerPatchyParticles->getPatchesInteractorInfo(inte.first).entryType    != "NonBondedPatches" or
                               (polymerPatchyParticles->getPatchesInteractorInfo(inte.first).entrySubType != "HelixExponential" and
                                polymerPatchyParticles->getPatchesInteractorInfo(inte.first).entrySubType != "HelixCosine")){
                                    System::log<System::CRITICAL>("[PatchPolymers] The patches of the DynamicallyBondedPatchyParticles"
                                                                  " must be NonBondedPatches of type HeliExponential or HelixCosine.");
                            }
                    }
                                }

                                if(isSurface){
                                        //Set type ids
                                        auto types = surfacePatchyParticles->getPatchesGlobalData()->getTypes();

                                        linkerType = types->getTypeId(linkerTypeName);

                                        //Check surface has the correct type
                                        for(auto inte : surfacePatchyParticles->getPatchesInteractors()){
                                                if(surfacePatchyParticles->getPatchesInteractorInfo(inte.first).entryType != "SurfacePatches"){
                                                        System::log<System::CRITICAL>("[PatchPolymers] The patches of the PatchyParticles must be SurfacePatches.");
                                                } else {
                                                        linker = std::dynamic_pointer_cast<typename Interactor::SingleInteractor<Potentials::SurfacePatches::Linker>>(inte.second)->getPotential();
                                                }
                                        }
                                }

                                //Fill the polymer group and the surface group

                                //Polymer group
                                {

                                        //Get the ids of the particles in the group
                                        auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                        access::location::cpu);

                                        auto polPPd = polymerPatchyParticles->getPatchesParticleData();

                                        auto patchesId = polPPd->getId(access::location::cpu,access::mode::read);
                                        auto parentId  = polPPd->getModelId(access::location::cpu,access::mode::read);

                                        std::vector<int> polymerIds;
                                        for(int i=0; i <this->pg->getNumberParticles(); i++){
                                                for(int j=0;j<polPPd->getNumParticles();j++){
                                                        if(pgIds[i] == parentId[j]){
                                                                polymerIds.push_back(patchesId[j]);
                                                        }
                                                }
                                        }

                                        //Create a new group with the polymer ids
                                        GroupUtils::selectors::ids selector(polymerIds);
                                        polymerGroup = std::make_shared<ParticleGroup>(selector,
                                                        polPPd,
                                                        "polymerGroup_"+pg->getName());

                                }

                                //If surface added create surfaceGroup
                                if(isSurface){

                                        //Get the ids of the particles in the group
                                        auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                        access::location::cpu);

                                        auto surfPPd = surfacePatchyParticles->getPatchesParticleData();

                                        auto patchesId = surfPPd->getId(access::location::cpu,access::mode::read);
                                        auto parentId  = surfPPd->getModelId(access::location::cpu,access::mode::read);

                                        std::vector<int> surfaceIds;
                                        for(int i=0; i<this->pg->getNumberParticles(); i++){
                                                for(int j=0;j<surfPPd->getNumParticles();j++){
                                                        if(pgIds[i] == parentId[j]){
                                                                surfaceIds.push_back(patchesId[j]);
                                                        }
                                                }
                                        }

                                        //Create a new group with the surface ids
                                        GroupUtils::selectors::ids selector(surfaceIds);
                                        surfaceGroup = std::make_shared<ParticleGroup>(selector,
                                                        surfacePatchyParticles->getPatchesParticleData(),
                                                        "surfaceGroup_"+pg->getName());
                                }

                                //Set up buffer
                                if(!isBufferSet){

                                        thrust::host_vector<int2> monomer2patches_h;
                                        thrust::host_vector<int>  patch2monomer_h;

                                        auto polPPd = polymerPatchyParticles->getPatchesParticleData();

                                        monomer2patches_h.resize(this->pd->getNumParticles(),{-1,-1});
                                        patch2monomer_h.resize(polPPd->getNumParticles(),-1);

                                        auto parentId   = polPPd->getModelId(access::location::cpu,access::mode::read);

                                        auto patchesId  = polPPd->getId(access::location::cpu,access::mode::read);
                                        auto patchesPos = polPPd->getPos(access::location::cpu,access::mode::read);

                                        for(int index=0; index<polPPd->getNumParticles(); index++){

                                                int currentParentId = parentId[index];
                                                int tpy = int(patchesPos[index].w);

                                                if(tpy == startType){
                                                        monomer2patches_h[currentParentId].y = patchesId[index];
                                                        patch2monomer_h[patchesId[index]]    = currentParentId;
                                                } else if(tpy == endType){
                                                        monomer2patches_h[currentParentId].x = patchesId[index];
                                                        patch2monomer_h[patchesId[index]]    = currentParentId;
                                                } else {
                                                        System::log<System::CRITICAL>("[PatchPolymers] The type of the patch (%d)"
                                                                        "is not the start (%d) or the end (%d) type.",
                                                                        tpy, startType, endType);
                                                }

                                        }

                                        for(int pId = 0; pId < this->pd->getNumParticles(); pId++){
                                                if(monomer2patches_h[pId].x == -1 || monomer2patches_h[pId].y == -1){
                                                        System::log<System::CRITICAL>("[PatchPolymers] The monomer with id %d"
                                                                        "does not have a start or end patch.", pId);
                                                }
                                        }

                                        for(int pId = 0; pId < polPPd->getNumParticles(); pId++){
                                                if(patch2monomer_h[pId] == -1){
                                                        System::log<System::CRITICAL>("[PatchPolymers] The patch with id %d"
                                                                        "does not have a monomer.", pId);
                                                }
                                        }

                                        patch2monomer   = patch2monomer_h;
                                        monomer2patches = monomer2patches_h;

                                        if(isSurface){

                                                thrust::host_vector<int> monomer2surfPatches_h;

                                                auto surfPPd = surfacePatchyParticles->getPatchesParticleData();

                                                monomer2surfPatches_h.resize(this->pd->getNumParticles(),-1);

                                                auto parentId   = surfPPd->getModelId(access::location::cpu,access::mode::read);

                                                auto patchesId  = surfPPd->getId(access::location::cpu,access::mode::read);
                                                auto patchesPos = surfPPd->getPos(access::location::cpu,access::mode::read);

                                                for(int index=0; index<surfPPd->getNumParticles(); index++){

                                                        int currentParentId = parentId[index];
                                                        int tpy = int(patchesPos[index].w);

                                                        if(tpy == linkerType){
                                                                monomer2surfPatches_h[currentParentId] = patchesId[index];
                                                        } else {
                                                                System::log<System::CRITICAL>("[PatchPolymers] The type of the patch (%d)"
                                                                                "is not linker type (%d).",
                                                                                tpy, linkerType);
                                                        }

                                                }

                                                for(int pId = 0; pId < this->pd->getNumParticles(); pId++){
                                                        if(monomer2surfPatches_h[pId] == -1){
                                                                System::log<System::CRITICAL>("[PatchPolymers] The monomer with id %d"
                                                                                "does not have a linker patch.", pId);
                                                        }
                                                }

                                                monomer2surfPatches = monomer2surfPatches_h;
                                        }

                                        //Set up buffer
                                        cudaError_t status = cudaMallocHost((void**)&maxEventsCPU, sizeof(int));
                                        if (status != cudaSuccess){
                                                System::log<System::CRITICAL>("[PatchPolymers] Error allocating memory on the CPU.");
                                        }
                                        maxEventsGPU.resize(1);

                                        eventBufferGPU.resize(this->pd->getNumParticles()*(bufferSize+6)); //6 is the number of events per particle
                                        lastBufferStep = this->gd->getFundamental()->getCurrentStep();

                                        isBufferSet = true;
                                }

                        }

                        ~PatchPolymers(){
                                System::log<System::DEBUG>("[PatchPolymers] Destroying patch polymers.");
                                if(isBufferSet){

                                        cudaDeviceSynchronize();
                                        resetBuffer(0);
                                        cudaDeviceSynchronize();

                                        cudaError_t status = cudaFreeHost(maxEventsCPU);
                                        if (status != cudaSuccess){
                                                System::log<System::CRITICAL>("[PatchPolymers] Error freeing memory on the CPU.");
                                        }

                                        monomer2patches.clear();
                                        monomer2patches.shrink_to_fit();

                                        patch2monomer.clear();
                                        patch2monomer.shrink_to_fit();

                                        monomer2surfPatches.clear();
                                        monomer2surfPatches.shrink_to_fit();

                                        nEventsPerMonomerGPU.clear();
                                        nEventsPerMonomerGPU.shrink_to_fit();

                                        eventBufferGPU.clear();
                                        eventBufferGPU.shrink_to_fit();

                                        maxEventsGPU.clear();
                                        maxEventsGPU.shrink_to_fit();

                                        previousState.clear();
                                        previousState.shrink_to_fit();

                                        previousSurfaceEnergy.clear();
                                        previousSurfaceEnergy.shrink_to_fit();

                                        isBufferSet = false;
                                }

                                System::log<System::DEBUG>("[PatchPolymers] Writing remaining events to file.");
                                writeEvents();
                        }

                        void applyStep(ullint step, cudaStream_t st) override {

                                if(!isBufferSet){
                                        System::log<System::CRITICAL>("[PatchPolymers] The buffer is not set.");
                                }

                                //Update buffer
                                if(step > lastBufferStep or firstStep){

                                        if(firstStep){

                                                auto polPPd = polymerPatchyParticles->getPatchesParticleData();

                                                thrust::host_vector<int4> previousState_h(polPPd->getNumParticles(),{-1,-1,-1,-1});
                                                previousState         = previousState_h;

                                                if(isSurface){
                                                        thrust::host_vector<real> previousSurfaceEnergy_h(this->pd->getNumParticles(),0.0);
                                                        previousSurfaceEnergy = previousSurfaceEnergy_h;
                                                }

                                                nEventsPerMonomerCPU.resize(this->pd->getNumParticles(),0);
                                                nEventsPerMonomerGPU = nEventsPerMonomerCPU;

                                                maxEventsCPU[0] = 0;
                                                cudaMemcpyAsync(thrust::raw_pointer_cast(maxEventsGPU.data()), maxEventsCPU,
                                                                sizeof(int), cudaMemcpyHostToDevice,st);

                                                cudaStreamSynchronize(st);

                                                if(isSurface){
                                                        surfacePatchyParticles->updatePatchyParticles(st);
                                                }
                                                cudaStreamSynchronize(st);

                                                firstStep = false;
                                        }

                                        auto polPPd = polymerPatchyParticles->getPatchesParticleData();

                                        int Nmonomers = this->pd->getNumParticles();
                                        int Npatches  = polPPd->getNumParticles();

                                        auto id2indexMonomers = this->pd->getIdOrderedIndices(access::location::gpu);
                                        auto id2indexPatches  = polPPd->getIdOrderedIndices(access::location::gpu);

                                        auto currentState     = polPPd->getState(access::location::gpu,access::mode::read).raw();

                                        int Nthreads = THREADS_PER_BLOCK;
                                        int Nblocks  = Nmonomers/Nthreads + ((Nmonomers%Nthreads)?1:0);

                                        if(!isSurface){
                                                PatchPolymers_ns::recordEvents
                                                        <<<Nblocks,Nthreads,0, st>>>(step,
                                                                        Nmonomers,
                                                                        bufferSize+6,
                                                                        id2indexMonomers,
                                                                        id2indexPatches,
                                                                        thrust::raw_pointer_cast(monomer2patches.data()),
                                                                        thrust::raw_pointer_cast(patch2monomer.data()),
                                                                        thrust::raw_pointer_cast(previousState.data()),
                                                                        currentState,
                                                                        thrust::raw_pointer_cast(eventBufferGPU.data()),
                                                                        thrust::raw_pointer_cast(nEventsPerMonomerGPU.data()),
                                                                        thrust::raw_pointer_cast(maxEventsGPU.data()));
                                        } else {

                                                auto surfPPd = surfacePatchyParticles->getPatchesParticleData();
                                                auto id2indexSurfPatches = surfPPd->getIdOrderedIndices(access::location::gpu);

                                                PatchPolymers_ns::recordEventsSurf
                                                        <<<Nblocks,Nthreads,0, st>>>(step,
                                                                        Nmonomers,
                                                                        bufferSize+6,
                                                                        id2indexMonomers,
                                                                        id2indexPatches,
                                                                        id2indexSurfPatches,
                                                                        thrust::raw_pointer_cast(monomer2patches.data()),
                                                                        thrust::raw_pointer_cast(patch2monomer.data()),
                                                                        thrust::raw_pointer_cast(monomer2surfPatches.data()),
                                                                        thrust::raw_pointer_cast(previousState.data()),
                                                                        currentState,
                                                                        thrust::raw_pointer_cast(previousSurfaceEnergy.data()),
                                                                        linker->getEnergyTransverser(),
                                                                        linker->getComputationalData(Computables(),st),
                                                                        surfaceEnergyThreshold,
                                                                        thrust::raw_pointer_cast(eventBufferGPU.data()),
                                                                        thrust::raw_pointer_cast(nEventsPerMonomerGPU.data()),
                                                                        thrust::raw_pointer_cast(maxEventsGPU.data()));

                                        }

                                        cudaMemcpyAsync(maxEventsCPU,
                                                        thrust::raw_pointer_cast(maxEventsGPU.data()),
                                                        sizeof(int), cudaMemcpyDeviceToHost,st);

                                        cudaStreamSynchronize(st);
                                        CudaCheckError();

                                        int maxEvents = maxEventsCPU[0];
                                        if(maxEvents > bufferSize){
                                                if(lastUpdateStep == std::numeric_limits<ullint>::max()){
                                                        System::log<System::DEBUG>("[PatchPolymers] Buffer full. "
                                                                        "Copying to CPU and resetting. "
                                                                        "Steps since last update: %llu",
                                                                        step);
                                                } else {
                                                        System::log<System::DEBUG>("[PatchPolymers] Buffer full. "
                                                                        "Copying to CPU and resetting. "
                                                                        "Steps since last update: %llu",
                                                                        step-lastUpdateStep);
                                                }

                                                resetBuffer(st);

                                                lastUpdateStep = step;
                                        }

                                        lastBufferStep = step;

                                }

                                if(step == lastUpdateStep){
                                        writeEvents();
                                }

                        }

        };

        //Static member initialization
        bool PatchPolymers::areParametersSet = false;
        bool PatchPolymers::isBufferSet = false;

        std::string PatchPolymers::startTypeName;
        std::string PatchPolymers::endTypeName;
        std::string PatchPolymers::linkerTypeName;

        int  PatchPolymers::startType;
        int  PatchPolymers::endType;
        int  PatchPolymers::linkerType;

        real PatchPolymers::surfaceEnergyThreshold;

        thrust::device_vector<int2> PatchPolymers::monomer2patches;
        thrust::device_vector<int>  PatchPolymers::patch2monomer;

        thrust::device_vector<int>  PatchPolymers::monomer2surfPatches;

        ullint PatchPolymers::bufferSize     = 0;

        ullint PatchPolymers::lastBufferStep = std::numeric_limits<ullint>::max();
        ullint PatchPolymers::lastUpdateStep = std::numeric_limits<ullint>::max();

        bool PatchPolymers::firstStep = true;

        thrust::host_vector<int>   PatchPolymers::nEventsPerMonomerCPU;
        thrust::device_vector<int> PatchPolymers::nEventsPerMonomerGPU;

        thrust::host_vector<PatchPolymers_ns::event>   PatchPolymers::eventBufferCPU;
        thrust::device_vector<PatchPolymers_ns::event> PatchPolymers::eventBufferGPU;

        int*                        PatchPolymers::maxEventsCPU;
        thrust::device_vector<int>  PatchPolymers::maxEventsGPU;

        thrust::device_vector<int4> PatchPolymers::previousState;
        thrust::device_vector<real> PatchPolymers::previousSurfaceEnergy;

        thrust::host_vector<int>                       PatchPolymers::nEventsPerMonomerRead;
        thrust::host_vector<PatchPolymers_ns::event>   PatchPolymers::PatchPolymers::eventBufferRead;


}}}}

REGISTER_SIMULATION_STEP(
    TopologicalMeasures,PatchPolymers,
    uammd::structured::SimulationStep::SimulationMeasures::PatchPolymers
)
