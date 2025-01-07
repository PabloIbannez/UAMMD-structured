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
        //Pp: polymerization positive direction
        //Pn: polymerization negative direction
        //Dp: depolymerization positive direction
        //Dn: depolymerization negative direction
        //S:  from bulk to surface
        //B:  from surface to bulk

        //Negative direction events (Pn and Dn) are ignored.
        //The reason is that the polymerization and depolymerization in the negative direction
        //can be inferred from the positive direction events.


        #pragma pack(push, 1) //Remove padding bytes
        struct event{
                ullint    step; //step in which the event occurs
                int       id;   //id of the particle (monomer)
                eventType type; //type of event
                int       info; //info about the event
        };
        #pragma pack(pop)

        __device__ inline void recordPolymerizationEvents(int id,
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

                //Note each thread is responsible for a monomer (not a patch)
                //Each monomer has two patches, one at the beginning and one at the end

                const int endPatchId = monomer2patches[id].x;
                const int begPatchId = monomer2patches[id].y;

                const int currEndState = currentState[id2indexPatches[endPatchId]].x;
                const int currBegState = currentState[id2indexPatches[begPatchId]].x;

                //Check endPatch
                // We can see the same event from the endPatch or the begPatch
                // I mean, a polymerization/depolymerization event in the endPatch appears as a polymerization/depolymerization event in the begPatch
                // The only difference is the info field, which is the patch id. But if we consider both options the same event is detected twice.
                // This is the reason why we only check the begPatch and all this code is commented.

                // We ensure each event is detected only once by checking the previous state of the begPatch

                //const int prevEndState = previousState[endPatchId].x;

                //if(prevEndState != currEndState){
                //	//Event detected
                //	if(prevEndState == int(-1)){
                //		//Polymerization
                //		event newEvent;

                //		newEvent.step = step;
                //		newEvent.id   = id;
                //		newEvent.type = Pn;
                //		newEvent.info = patch2monomer[currEndState];

                //		UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                //							   bufferSize, Nmonomers,
                //							   nEventsPerMonomer[id], id,
                //							   newEvent);

                //		nEventsPerMonomer[id]++;
                //	}

                //	if(currEndState == int(-1)){
                //		//Depolymerization
                //		event newEvent;

                //		newEvent.step = step;
                //		newEvent.id   = id;
                //		newEvent.type = Dn;
                //		newEvent.info = patch2monomer[prevEndState];

                //		UAMMD_SET_2D_ROW_MAJOR(eventBuffer,
                //							   bufferSize, Nmonomers,
                //							   nEventsPerMonomer[id], id,
                //							   newEvent);

                //		nEventsPerMonomer[id]++;
                //	}
                //}

                //Check begPatch
                const int prevBegState = previousState[begPatchId].x;

                if(prevBegState != currBegState){
                        //Event detected

                        if(prevBegState == int(-1)){
                                //Polymerization
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
                                //Depolymerization
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

                recordPolymerizationEvents(id,
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

                recordPolymerizationEvents(id,
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

            std::ofstream outputFileEvents;
            std::ofstream outputFilePolymerDistribution;
            std::ofstream outputFileSurfacePolymerDistribution;
            std::ofstream outputFileFractionMonomersInSurface;
            std::ofstream outputFileTransitionMatrix;

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

            static bool writeEvents;

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

            //Measuring

            int nMonomers;
            int idOffset;

            int maxLen;

	        std::vector<int>  monomer2polymer;
	        std::vector<bool> isMonomerSurface;

            int nextPolymerId = 0;
	        std::map<int,std::vector<int>> polymers;
	        std::set<int>              ringPolymers;

	        //struct birthDeath{
	        //	ullint birth;
	        //	ullint death;
	        //};

	        //std::vector<birthDeath>   polymersLife;
	        std::set<int>             polymersAlive;

	        std::map<std::pair<std::string,std::string>,ullint> transitions;

	        struct stateType{
	        	ullint step;
	        	bool surface;
	        	bool ring;
	        	int  left;
	        	int  right;
	        };

            std::vector<stateType> monomerState;

	        stateType getMonomerState(ullint step, int id){

	        	stateType state;
	        	state.step = step;

	        	int polymerId = monomer2polymer[id];
	        	if(polymerId == -1){

	        		state.left  = 0;
	        		state.right = 0;

	        		if(isMonomerSurface[id]){
	        			state.surface = true;
	        		} else {
	        			state.surface = false;
	        		}

	        	} else {

	        		const std::vector<int>& polymer = polymers[polymerId];

	        		//Find the position of the monomer in the polymer
	        		int index = -1;
	        		for(int i=0; i<polymer.size(); i++){
	        			if(polymer[i] == id){
	        				index = i;
	        				break;
	        			}
	        		}

	        		if(index == -1){
                        System::log<System::CRITICAL>("[PatchPolymers] Monomer not found in the polymer");
	        		}

	        		//Count the number of monomers left
	        		int left = 0;
	        		for(int i=0; i<index; i++){
	        				left++;
	        		}
	        		int right = 0;
	        		for(int i=index+1; i<polymer.size(); i++){
	        				right++;
	        		}

	        		state.left  = left;
	        		state.right = right;

	        		//If some monomer is on the surface, the polymer is on the surface
	        		state.surface = false;
	        		for(int i=0; i<polymer.size(); i++){
	        			if(isMonomerSurface[polymer[i]]){
	        				state.surface = true;
	        				break;
	        			}
	        		}

	        	}

	        	if(isRing(id)){
	        		state.ring = true;
	        	} else {
	        		state.ring = false;
	        	}

	        	return state;
	        }

		    bool isPolymerSurface(int pol){
		    	std::vector<int> polymer = polymers[pol];
		    	bool surf = false;
		    	for(int i=0; i<polymer.size(); i++){
		    		if(isMonomerSurface[polymer[i]]){
		    			surf = true;
		    			break;
		    		}
		    	}

		    	return surf;
		    }

	        std::string state2key(const stateType& state){
	        	std::string key = "";

	        	int left  = state.left;
	        	int right = state.right;

	        	int len = right+left;

	        	bool ring    = state.ring;
	        	bool surface = state.surface;

	        	if(surface){
	        		key += "S_"; //Surface
	        	} else {
	        		key += "B_"; //Bulk
	        	}

	        	if(len <= maxLen){
	        		if(ring){
	        			key += "R_"; //Ring
	        			key = key + std::to_string(right + left);
	        		} else {
	        			key += "L_"; //Linear
	        			key = key + std::to_string(left) + "_" + std::to_string(right);
	        		}
	        	} else {
	        		if(ring){
	        			key += "RM"; //Ring
	        		} else {
	        			key += "LM"; //Linear
	        		}
	        	}

	        	return key;
	        }

            //Measuring functions

	        void processPolymerization(const PatchPolymers_ns::event& e){

	        	ullint step = e.step;

	        	int resultingPolymerId; //Id of the resulting polymer

	        	int id1 = e.id;
	        	int id2 = e.info;

	        	//id1 --> id2, order matters !!!

	        	int polymer1 = monomer2polymer[id1];
	        	int polymer2 = monomer2polymer[id2];

	        	if(polymer1 == -1 and polymer2 == -1){
	        		//Two monomers to a new polymer
	        		resultingPolymerId = createPolymer(step,id1,id2);
	        	} else {
	        		if(polymer1 == polymer2){
	        			//Polymer merges itself, a ring is formed
	        			resultingPolymerId = createRing(id1,id2);
	        		} else {
	        			if        (polymer1 == -1){
	        				//Merge monomer 1 to polymer 2
	        				resultingPolymerId = mergeMonomerPolymer(step,id1,id2);
	        			} else if (polymer2 == -1){
	        				//Merge polymer 1 to monomer 2
	        				resultingPolymerId = mergePolymerMonomer(id1,id2);
	        			} else {
	        				//Two polymers merge
	        				resultingPolymerId = mergePolymers(step,id1,id2);
	        			}
	        		}
	        	}
	        }

	        void processDepolymerization(const PatchPolymers_ns::event& e){

	        	ullint step = e.step;

	        	int resultingPolymerId1; //Id of the resulting polymer
	        	int resultingPolymerId2; //Id of the resulting polymer

	        	int id1 = e.id;
	        	int id2 = e.info;

	        	//id1 -/-> id2, order matters !!!

	        	int polymer1 = monomer2polymer[id1];
	        	int polymer2 = monomer2polymer[id2];

	        	if(polymer1 == -1 and polymer2 == -1){
	        		//This should not happen
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to depolymerize a monomer from a monomer");
	        	} else {
	        		if(isRing(id1)){
	        			//Polymer splits itself, a ring is broken
	        			int2 resultingPolymers = destroyRing(id1,id2);
	        			resultingPolymerId1 = resultingPolymers.x;
	        			resultingPolymerId2 = resultingPolymers.y;
	        		} else {
	        			if        (polymer1 == -1){
	        				//This should not happen
                            System::log<System::CRITICAL>("[PatchPolymers] Trying to depolymerize a monomer from a polymer");
	        			} else if (polymer2 == -1){
	        				//This should not happen
                            System::log<System::CRITICAL>("[PatchPolymers] Trying to depolymerize a polymer from a monomer");
	        			} else {
	        				//Two polymers split
	        				int2 resultingPolymers = splitPolymer(step,id1,id2);
	        				resultingPolymerId1 = resultingPolymers.x;
	        				resultingPolymerId2 = resultingPolymers.y;
	        			}
	        		}
	        	}
	        }

	        void processBulkToSurface(const PatchPolymers_ns::event& e){

	        	ullint step = e.step;
	        	int id = e.id;

	        	if(isMonomerSurface[id]){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to move a monomer to the surface when it is already there");
	        	} else {
	        		isMonomerSurface[id] = true;
	        	}
	        }

	        void processSurfaceToBulk(const PatchPolymers_ns::event& e){

	        	ullint step = e.step;
	        	int id = e.id;

	        	if(!isMonomerSurface[id]){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to move a monomer to the bulk when it is already there");
	        	} else {
	        		isMonomerSurface[id] = false;
	        	}
	        }

            //Polymerization functions

	        int createPolymer(ullint step, int id1,int id2){
	        	//Polymer is created only when two monomers are connected
	        	//Check both monomers don't belong to a polymer
	        	if(monomer2polymer[id1] != -1 or
	        	   monomer2polymer[id2] != -1){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a polymer with a monomer that already belongs to a polymer");
	        	}

	        	////////////////////////////////
	        	int polymerId = nextPolymerId;
                nextPolymerId++;
	        	//Create an empty polymer
                //Check if polymerId is already in use
                if(polymers.find(polymerId) != polymers.end()){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a polymer with an id that is already in use");
                }
                polymers[polymerId] = {id1,id2};

	        	monomer2polymer[id1] = polymerId;
	        	monomer2polymer[id2] = polymerId;

	        	//polymersLife.push_back({step,0});
                //Check if polymerId is already in use
                if(polymersAlive.find(polymerId) != polymersAlive.end()){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to add a polymer to the alive polymers that is already there");
                }
	        	polymersAlive.insert(polymerId);
	        	//////////////////////////////////

	        	return polymerId;
	        }

	        int createPolymer(ullint step, std::vector<int> ids){

	        	//Check all monomers don't belong to a polymer
	        	for(int id:ids){
	        		if(monomer2polymer[id] != -1){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to create a polymer with a monomer that already belongs to a polymer");
	        		}
	        	}

	        	//Check that size of ids is at least 2
	        	if(ids.size() < 2){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a polymer with less than two monomers");
	        	}

	        	//Create polymer
	        	int polymerId = createPolymer(step,ids[0],ids[1]);
	        	std::vector<int>& polymer = polymers[polymerId];

	        	//Add the rest of monomers to the polymer
	        	for(int i=2;i<ids.size();i++){
	        		int id = ids[i];
	        		polymer.push_back(id);
	        	}

	        	//Update monomer2polymer
	        	for(int id:polymer){
	        		monomer2polymer[id] = polymerId;
	        	}

	        	return polymerId;
	        }

	        std::vector<int> destroyPolymer(ullint step, int id1){
	        	int polymerId = monomer2polymer[id1];

	        	//If polymer is a ring, error
	        	if(ringPolymers.find(polymerId) != ringPolymers.end()){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to destroy a ring polymer");
	        	}

	        	std::vector<int> polymer = polymers[polymerId];

	        	////////////////////////////////
	        	for(int i=0; i<polymer.size(); i++){
	        		monomer2polymer[polymer[i]] = -1;
	        	}
	        	//polymersLife[polymerId].death = step;
	        	polymersAlive.erase(polymerId);
                polymers.erase(polymerId);
	        	////////////////////////////////

	        	return polymer;
	        }

	        int createRing(int id1,int id2){
	        	//Ring is created only when one polymer is connected to itself
	        	//Check both monomers belong to the same polymer

	        	if(monomer2polymer[id1] != monomer2polymer[id2]){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a ring with monomers that belong to different polymers");
	        	}

	        	//Check that id1 is the last monomer in the polymer and id2 is the first
	        	if(polymers[monomer2polymer[id1]].back() != id1 or
	        	   polymers[monomer2polymer[id2]].front() != id2){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a ring with monomers that are not the last and first in the polymer");
	        	}

	        	//Check if polymer is already a ring
	        	if(ringPolymers.find(monomer2polymer[id1]) != ringPolymers.end()){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to create a ring with a polymer that is already a ring");
	        	}

	        	int polymerId = monomer2polymer[id1];
	        	ringPolymers.insert(polymerId);

	        	return polymerId;
	        }

	        int2 destroyRing(int id1,int id2){
	        	//Ring is destroyed only when one polymer is connected to itself
	        	//Check both monomers belong to the same polymer

	        	//Check that both are ring
	        	if(not isRing(id1) or not isRing(id2)){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to destroy a ring with monomers that don't belong to a ring");
	        	}

	        	if(monomer2polymer[id1] != monomer2polymer[id2]){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to destroy a ring with monomers that belong to different polymers");
	        	}

	        	int polymerId = monomer2polymer[id1];
                ringPolymers.erase(polymerId);

	        	//Update polymer so id1 is the last monomer and id2 is the first
	        	std::vector<int> polymer = polymers[polymerId];

	        	int index1 = std::find(polymer.begin(), polymer.end(), id1) - polymer.begin();
	        	int index2 = std::find(polymer.begin(), polymer.end(), id2) - polymer.begin();

	        	//Check that indices are consecutive or the last and first
	        	bool consecutive  = (index1 == index2+1) or (index1 == index2-1);
	        	bool lastAndFirst = (index1 == polymer.size()-1) and (index2 == 0);

	        	if(not (consecutive or lastAndFirst)){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to destroy a ring with monomers that are not consecutive or the last and first");
	        	}

	        	std::rotate(polymer.begin(), std::find(polymer.begin(), polymer.end(), id2), polymer.end());

	        	////Print, for checking
	        	//std::cout << "Ring destroyed in bond " << id1 << "-" << id2 << std::endl;
	        	//std::cout << "Old polymer: ";
	        	//for(int i=0; i<polymers[polymerId].size(); i++){
	        	//	std::cout << polymers[polymerId][i] << " ";
	        	//}
	        	//std::cout << std::endl;
	        	//std::cout << "New polymer: ";
	        	//for(int i=0; i<polymer.size(); i++){
	        	//	std::cout << polymer[i] << " ";
	        	//}
	        	//std::cout << std::endl;
	        	//std::cin.get();

	        	polymers[polymerId] = polymer;

	        	return {polymerId,polymerId};
	        }

	        bool isRing(int id){
	        	int polymerId = monomer2polymer[id];
	        	if(ringPolymers.find(polymerId) != ringPolymers.end()){
	        		return true;
	        	}
	        	return false;
	        }

	        int mergePolymers(ullint step,int id1, int id2){
	        	int polymerId1 = monomer2polymer[id1];
	        	int polymerId2 = monomer2polymer[id2];

	        	//Check both monomers belong to different polymers
	        	if(polymerId1 == -1 or
	        	   polymerId2 == -1){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to merge polymers with monomers that don't belong to a polymer");
	        	}

	        	if(polymerId1 == polymerId2){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to merge a polymer with itself");
	        	}

	        	//Check that id1 is the last monomer in the polymer and id2 is the first
	        	if(polymers[polymerId1].back()  != id1 or
	        	   polymers[polymerId2].front() != id2){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to merge polymers with monomers that are not the last and first in the polymer");
	        	}

	        	//Destroy polymer 2
	        	std::vector<int> polymer2 = destroyPolymer(step,id2);

	        	//Add polymer 2 to polymer 1
	        	polymers[polymerId1].insert(polymers[polymerId1].end(),polymer2.begin(),polymer2.end());

	        	//Update monomer2polymer
	        	for(int i=0; i<polymer2.size(); i++){
	        		monomer2polymer[polymer2[i]] = polymerId1;
	        	}

	        	return polymerId1;
	        }

	        int mergePolymerMonomer(int id1, int id2){ //Note step is not needed
	        	int polymerId1 = monomer2polymer[id1];
	        	int polymerId2 = monomer2polymer[id2];

	        	//Check id1 is a polymer and id2 is a monomer
	        	if(polymerId1 == -1 or
	        	   polymerId2 != -1){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to merge a polymer with a monomer that don't belong to a polymer");
	        	}

	        	//Add monomer to polymer
	        	polymers[polymerId1].push_back(id2);
	        	monomer2polymer[id2] = polymerId1;

	        	return polymerId1;
	        }

	        int mergeMonomerPolymer(ullint step, int id1, int id2){
	        	int polymerId1 = monomer2polymer[id1];
	        	int polymerId2 = monomer2polymer[id2];

	        	//Check id1 is a monomer and id2 is a polymer
	        	if(polymerId1 != -1 or
	        	   polymerId2 == -1){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to merge a monomer with a polymer that don't belong to a polymer");
	        	}

	        	//Destroy polymer 2
	        	std::vector<int> polymerNew = destroyPolymer(step,id2);
	        	//Add monomer to polymer, at the beginning
	        	polymerNew.insert(polymerNew.begin(),id1);

	        	int polymerId = createPolymer(step,polymerNew);

	        	return polymerId;
	        }

	        int2 splitPolymer(ullint step, int id1, int id2){
	        	int polymerId1 = monomer2polymer[id1];
	        	int polymerId2 = monomer2polymer[id2];

	        	//Check both are polymers
	        	if(polymerId1 == -1 or
	        	   polymerId2 == -1){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers with monomers that don't belong to a polymer");
	        	}

	        	//Check both belong to the same polymer
	        	if(polymerId1 != polymerId2){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers with monomers that belong to different polymers");
	        	}

	        	std::vector<int> polymer = polymers[polymerId1];

	        	//Find the index of the monomers
	        	int index1 = -1;
	        	int index2 = -1;
	        	for(int i=0; i<polymer.size(); i++){
	        		if(polymer[i] == id1){
	        			index1 = i;
	        		}
	        		if(polymer[i] == id2){
	        			index2 = i;
	        		}
	        	}

	        	//Check both are consecutive
	        	if(index1 != (index2-1)){
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers with monomers that are not consecutive");
	        	}

	        	bool is1Monomer = (index1 == 0); //Will id1 be a monomer?
	        	bool is2Monomer = (index2 == (polymer.size()-1)); //Will id2 be a monomer?

	        	//Cases:
	        	// 1) 1 -> monomer, 2->monomer
	        	// 2) 1 -> monomer, 2->polymer
	        	// 3) 1 -> polymer, 2->monomer
	        	// 4) 1 -> polymer, 2->polymer

	        	int resultingPolymerId1; //Id of the resulting polymer
	        	int resultingPolymerId2; //Id of the resulting polymer
	        													 //
	        	if        (is1Monomer and is2Monomer){
	        		std::vector<int> tmp = destroyPolymer(step,id1);
	        		//Check that pol2 is a vector with two elements, id1 and id2
	        		if(tmp.size() != 2 or
	        			 tmp[0] != id1 or
	        			 tmp[1] != id2){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers (1)");
	        		}

	        		resultingPolymerId1 = -1;
	        		resultingPolymerId2 = -1;

	        	} else if (is1Monomer and not is2Monomer){
	        		//Check id1 is the first monomer of the polymer
	        		if(index1 != 0){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymer, but id1 is not the first monomer");
	        		}
	        		//Destroy polymer
	        		std::vector<int> tmp = destroyPolymer(step,id1);
	        		resultingPolymerId1  = -1;

	        		//Check the first element of pol2 is id1 and the second is id2
	        		if(tmp[0] != id1 or
	        		   tmp[1] != id2){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers (2)");
	        		}
	        		//Remove id1 from pol2
	        		tmp.erase(tmp.begin());

	        		//Create polymer, starting with id2
	        		resultingPolymerId2 = createPolymer(step,tmp);
	        	} else if (not is1Monomer and is2Monomer){
	        		//Check id2 is the last monomer of the polymer
	        		if(index2 != (polymer.size()-1)){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymer, but id2 is not the last monomer");
	        		}

	        		std::vector<int>& polymer = polymers[polymerId1];

	        		//Remove last element of pol1
	        		polymer.pop_back();

	        		//Update monomer2polymer
	        		monomer2polymer[id2] = -1;

	        		resultingPolymerId1 = polymerId1;
	        		resultingPolymerId2 = -1;
	        	} else if (not is1Monomer and not is2Monomer){

	        		std::vector<int>& polymer = polymers[polymerId1];
	        		std::vector<int> tmp;

	        		//Remove all elements from pol1 after id1
	        		for(int i=index1+1; i<polymer.size(); i++){
	        			tmp.push_back(polymer[i]);
	        			monomer2polymer[polymer[i]] = -1;
	        		}

	        		//Remove all elements from pol1 after id1
	        		polymer.erase(polymer.begin()+index1+1,polymer.end());

	        		//Check the last element of polymer is id1 and the first of tmp is id2
	        		if(polymer.back() != id1 or
	        		   tmp[0] != id2){
                        System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers (3)");
	        		}

	        		//Create polymer, starting with id2
	        		resultingPolymerId1 = polymerId1;
	        		resultingPolymerId2 = createPolymer(step,tmp);

	        	} else {
                    System::log<System::CRITICAL>("[PatchPolymers] Trying to split polymers (4)");
	        	}

	        	return {resultingPolymerId1,resultingPolymerId2};
	        }

            //Event processing functions

            void processEvents(const ullint& step){

                //Get the ids of the particles in the group
                //Here is the place where we particularize the group of particles
                auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                           access::location::cpu);

                eventWriteBuffer.resize(0);

                for(int i=0; i <this->pg->getNumberParticles(); i++){
                        int id = pgIds[i]; //Note that here the idOffset is not subtracted.
                                           //This is because the events are stored with the original ids are
                                           //shared between all PatchPolymers instances.
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

                // Write events to file
                if(writeEvents){
                        //Write the events to the file in binary format
                        for(auto e : eventWriteBuffer){
                                outputFileEvents.write((char*)&e, sizeof(PatchPolymers_ns::event));
                        }
                }

                ///////////////////////////
                // Start events analysis //
                ///////////////////////////

                for(auto e : eventWriteBuffer){
                    int type = e.type;
                    e.id -= idOffset;

                    System::log<System::DEBUG>("[PatchPolymers] Processing event, "
                                                 "step: %llu, id: %d, type: %d (0: Pp, 1: Pn, 2: Dp, 3: Dn, 4: S, 5: B), info: %d. Instance: %s",
                                                 e.step,e.id,type,e.info,this->name.c_str());

		            if       (type == PatchPolymers_ns::eventType::Pp){ //polymerization positive direction
                        System::log<System::DEBUG>("[PatchPolymers] Processing polymerization event.");
                        e.info -= idOffset; //The info is the second monomer
		            	processPolymerization(e);
		            } else if(type == PatchPolymers_ns::eventType::Dp){ //depolymerization positive direction
                        System::log<System::DEBUG>("[PatchPolymers] Processing depolymerization event.");
                        e.info -= idOffset; //The info is the second monomer
		            	processDepolymerization(e);
		            } else if(type == PatchPolymers_ns::eventType::S){ //from bulk to surface
                        System::log<System::DEBUG>("[PatchPolymers] Processing bulk to surface event.");
		            	processBulkToSurface(e);
		            } else if(type == PatchPolymers_ns::eventType::B){ //from surface to bulk
                        System::log<System::DEBUG>("[PatchPolymers] Processing surface to bulk event.");
		            	processSurfaceToBulk(e);
		            } else {
                        System::log<System::CRITICAL>("[PatchPolymers] Event type not recognized.");
		            }

                    //At this point, the event has been processed
                    stateType prevState = monomerState[e.id];
                    stateType newState  = getMonomerState(e.step,e.id);

				    ullint stepsInPrevState = newState.step - prevState.step;

                    auto ppKey = std::make_pair(state2key(prevState),state2key(prevState));
                    auto pnKey = std::make_pair(state2key(prevState),state2key(newState));

                    if (transitions.find(ppKey) == transitions.end()){
                        transitions[ppKey] = 0;
                    }

                    if (transitions.find(pnKey) == transitions.end()){
                        transitions[pnKey] = 0;
                    }

				    transitions[ppKey] += stepsInPrevState;
				    transitions[pnKey] += 1;

                    //Update the state of the monomer
                    monomerState[e.id] = newState;
                }
            }

            void resetBuffer(cudaStream_t st){
                // Reset the buffer copying the events to the CPU and resetting the GPU buffer
                // The array eventBufferRead is used to store the events in the CPU
                // The array nEventsPerMonomerRead is used to store the number of events per monomer in the CPU

                eventBufferRead       = eventBufferGPU;
                nEventsPerMonomerRead = nEventsPerMonomerGPU;

                /////////////////////////////////////////////

                //We reset the maxEvents to 0
                maxEventsCPU[0] = 0;
                cudaMemcpyAsync(thrust::raw_pointer_cast(maxEventsGPU.data()), maxEventsCPU,
                                sizeof(int), cudaMemcpyHostToDevice,st);

                //We reset the events counter per monomer to 0
                std::fill(nEventsPerMonomerCPU.begin(),nEventsPerMonomerCPU.end(),0);
                nEventsPerMonomerGPU = nEventsPerMonomerCPU;

                // WARNING: We do not reset the eventBufferGPU. We just overwrite the events in the buffer.
                // We get this behavior because we are resetting nEventsPerMonomerGPU to 0.
                // This array tells us where to write the next event in the buffer.
                // This is done for performance reasons. We do not need to reset the whole buffer, just the counter.

                cudaStreamSynchronize(st);

            }

            //Write functions

            void openFile(std::string name, std::ofstream& file){
                std::string filePath = this->outputFilePath + name;

                bool isFileEmpty = Backup::openFile(this->sys, filePath, file, true);

                if(not isFileEmpty){
                        System::log<System::MESSAGE>("[PatchPolymers] Output file (%s) is not empty. Appending data.", filePath.c_str());
                }
            }

            void writePolymerDistribution(const ullint& step){

			    int localMaxPolymerSize = 0;
			    for(std::set<int>::iterator it = polymersAlive.begin(); it != polymersAlive.end(); it++){
			    	int polyd = *it;
			    	int polSize = polymers[polyd].size()-1;
			    	localMaxPolymerSize = std::max(localMaxPolymerSize,polSize);
			    }

			    std::vector<int> distribution(localMaxPolymerSize+1,0);
			    std::vector<int> distributionSurface(localMaxPolymerSize+1,0);

			    //Count the number of monomers which are not in a polymer
			    for(int i=0; i<nMonomers; i++){
			    	if(monomer2polymer[i] == -1){
			    	if(not isMonomerSurface[i]){
			    		distribution[0]++;
			    	} else {
			    		distributionSurface[0]++;
			    	}
			    	}
			    }

			    //Iterate over polymers allive (set)
			    for(std::set<int>::iterator it = polymersAlive.begin(); it != polymersAlive.end(); it++){
			    	int polyd = *it;
			    	int polSize = polymers[polyd].size()-1;
			    	if(not isPolymerSurface(polyd)){
			    		distribution[polSize]++;
			    	} else {
			    		distributionSurface[polSize]++;
			    	}
			    }

                //Write the polymer distribution to the file outputFilePolymerDistribution
                outputFilePolymerDistribution << "# " << step << std::endl;
                for(int i=0; i<distribution.size(); i++){
                    outputFilePolymerDistribution << i << " " << distribution[i] << std::endl;
                }

                if(isSurface){
                    //Write the surface polymer distribution to the file outputFileSurfacePolymerDistribution
                    outputFileSurfacePolymerDistribution << "# " << step << std::endl;
                    for(int i=0; i<distributionSurface.size(); i++){
                        outputFileSurfacePolymerDistribution << i << " " << distributionSurface[i] << std::endl;
                    }
                }
            }

            void writeFractionMonomersInSurface(const ullint& step){
                if(not isSurface){
                    return;
                }

			    int nMonomersInSurface = 0;
			    for(int i=0; i<nMonomers; i++){
			    	if(isMonomerSurface[i]){
			    		nMonomersInSurface++;
			    	}
			    }

			    double fractionMonomersInSurface = (double)nMonomersInSurface/nMonomers;

                //Write the fraction of monomers in the surface to the file outputFileFractionMonomersInSurface
                outputFileFractionMonomersInSurface << step << " " << fractionMonomersInSurface << std::endl;
            }

            void writeTransitionMatrix(const ullint& step){
                //Write the transition matrix to the file outputFileTransitionMatrix
                outputFileTransitionMatrix << "# " << step << std::endl;
			    for(auto it=transitions.begin(); it!=transitions.end(); it++){
			    	std::string fromState = it->first.first;
			    	std::string toState   = it->first.second;

			    	ullint nEvents = it->second;

			    	outputFileTransitionMatrix << fromState << " " << toState << " " << nEvents << std::endl;
			    }
            }

        public:

            PatchPolymers(std::shared_ptr<ParticleGroup>             pg,
                          std::shared_ptr<IntegratorManager> integrator,
                          std::shared_ptr<ForceField>                ff,
                          DataEntry& data,
                          std::string name):SimulationStepBase(pg,integrator,ff,data,name){

                    //Read the output file path
                    outputFilePath = data.getParameter<std::string>("outputFilePath");

                    if(!areParametersSet){
                            startTypeName = data.getParameter<std::string>("startType","S");
                            endTypeName   = data.getParameter<std::string>("endType"  ,"E");
                            linkerTypeName= data.getParameter<std::string>("linkerType","L");

                            surfaceEnergyThreshold = data.getParameter<real>("surfaceEnergyThreshold",0.0);

                            writeEvents = data.getParameter<bool>("writeEvents",false);

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

                    //Not common parameters

                    maxLen = data.getParameter<int>("maxLen",25);

                    //

                    ullint newBufferSize = data.getParameter<ullint>("bufferSize");
                    if(newBufferSize != bufferSize and bufferSize != 0){
                            //Error
                            System::log<System::CRITICAL>("[PatchPolymers] Buffer step is not the same for all the instances of the class");
                    } else {
                            bufferSize = newBufferSize;
                    }

                    // Set up nMonomers and idOffset
                    auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                               access::location::cpu);

                    nMonomers = this->pg->getNumberParticles();

                    // Determine the idOffset
                    std::vector<int> ids(pgIds,pgIds+nMonomers);
                    std::sort(ids.begin(),ids.end()); //Sort the ids in ascending order

                    for(int i=1; i<nMonomers; i++){
                        if(ids[i] != ids[i-1]+1){
                            System::log<System::CRITICAL>("[PatchPolymers] The ids of the particles in the group are not consecutive.");
                        }
                    }

                    idOffset = ids[0];

                    //Init monomer2polymer and isMonomerSurface
                    monomer2polymer.resize(nMonomers,-1);
                    isMonomerSurface.resize(nMonomers,false);
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
                processEvents(0); //Write remaining events to file
            }

            //We have to override the init function.
            void init(cudaStream_t st) override{

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
                                   (polymerPatchyParticles->getPatchesInteractorInfo(inte.first).entrySubType != "COSRAP")){
                                        System::log<System::CRITICAL>("[PatchPolymers] The patches of the DynamicallyBondedPatchyParticles"
                                                                      " must be NonBondedPatches of type COSRAP.");
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
                                        linker = std::dynamic_pointer_cast<typename Interactor::SingleInteractor<Potentials::SurfacePatches::Linker>>
                                        (inte.second)->getPotential();
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

                    System::log<System::MESSAGE>("[PatchPolymers] Instance %s is setting up the buffer.",this->name.c_str());

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

                    // COMMENT ABOUT THE SIZE OF THE BUFFER
                    // You can see that the size of the buffer is bufferSize+6.
                    // 6 comes from the number of event types. Summing 6 to the bufferSize we
                    // ensure that the buffer is big enough to store all the events. Somehow it is assuming
                    // that a particle can have all the events at the same time. This is not true
                    // (I think the maximum would be 4, or 2 since we are ignoring the negative direction events).
                    // But in any case, the buffer is big enough to store all the events and we do not expect out of bounds errors.

                    eventBufferGPU.resize(this->pd->getNumParticles()*(bufferSize+6)); //6 is the number of event types
                    lastBufferStep = this->gd->getFundamental()->getCurrentStep();

                    isBufferSet = true;

                    System::log<System::MESSAGE>("[PatchPolymers] Buffer set up.");
                }

                //Open the output files

                System::log<System::MESSAGE>("[PatchPolymers] Opening output files.");

                if(writeEvents){
                    openFile("Events.dat",outputFileEvents);
                }

                openFile("PolymerDistribution.dat",outputFilePolymerDistribution);
                openFile("TransitionMatrix.dat",outputFileTransitionMatrix);

                if(isSurface){
                    openFile("SurfacePolymerDistribution.dat",outputFileSurfacePolymerDistribution);
                    openFile("FractionMonomersInSurface.dat",outputFileFractionMonomersInSurface);
                }

                //Init monomerState
                {
                    System::log<System::MESSAGE>("[PatchPolymers] Initializing monomer state.");

                    // Set up nMonomers and idOffset
                    monomerState.resize(nMonomers);

                    auto pgIds = this->pg->getPropertyIterator(this->pd->getId(access::location::cpu, access::mode::read).begin(),
                                                               access::location::cpu);

                    ullint currentStep = this->gd->getFundamental()->getCurrentStep();
                    for(int i=0; i<nMonomers; i++){
                        int id = pgIds[i]-idOffset;
                        stateType state = getMonomerState(currentStep,id);
                        monomerState[id] = state;
                    }
                }

                System::log<System::MESSAGE>("[PatchPolymers] Initialization done.");
            }

            void update(ullint step, cudaStream_t st) override {

                    //Update is called at each integration step. We record the events at each step (we record ALL the events).
                    //Events are stored in a buffer. When the buffer is full, we write the events to a file and reset the buffer.
                    //Events buffer is considered full when the maximum number of events is reached for a SINGLE monomer.

                    if(!isBufferSet){
                        System::log<System::CRITICAL>("[PatchPolymers] The buffer is not set.");
                    }

                    //Update buffer
                    if(step > lastBufferStep or firstStep){ // Comparing with lastBufferStep ensure that just one update is done across all the instances of the class.
                    // The different instances can work over different groups of particles. But this group is only taken into account when writing the events to the file.
                    // All the instances share the same buffer and the same buffer size, so we have to ensure that only one update is done across all the instances.

                        if(firstStep){

                            //Init stuff. Here we assume that at beginning there are no polymers or patches attached to the surface.

                            auto polPPd = polymerPatchyParticles->getPatchesParticleData();

                            thrust::host_vector<int4> previousState_h(polPPd->getNumParticles(),{-1,-1,-1,-1});
                            previousState         =   previousState_h;

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
                        System::log<System::DEBUG>("[PatchPolymers] Max number of events for a single monomer: %d", maxEvents);
                        //Here we check if we have reached the maximum number of events for a single monomer
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
                                // We can access the buffer in the CPU and write the events to the file.
                                // The information is stored in the CPU arrays eventBufferRead and nEventsPerMonomerRead.
                        }
                        lastBufferStep = step;

                    }

                    if(step == lastUpdateStep){
                            // This means that the buffer has been reset in this step.
                            // It is ensured that the buffer information is stored in the CPU arrays eventBufferRead and nEventsPerMonomerRead.
                            processEvents(step);
                    }

            }

            void applyStep(ullint step, cudaStream_t st) override {

                if(lastUpdateStep != step){
                    resetBuffer(st);
                    processEvents(step);
                    lastUpdateStep = step;

                    //This is SUBTLE. We only process the events in the instance that reset the buffer.
                    //The other instances will process the events in their update function.
                    //The execution scheme of the simulation steps is the following:
                    // 1: update
                    // 1: applyStep
                    // 2: update
                    // 2: applyStep
                    // 3: update
                    // 3: applyStep
                    // ...
                }


                //Write the data to the file
                writePolymerDistribution(step);
                writeFractionMonomersInSurface(step);
                writeTransitionMatrix(step);
            }

    };

    //Static member initialization
    bool PatchPolymers::areParametersSet = false;
    bool PatchPolymers::isBufferSet = false;

    bool PatchPolymers::writeEvents = false;

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
