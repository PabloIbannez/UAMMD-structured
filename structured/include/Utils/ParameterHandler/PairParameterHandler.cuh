#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Utils/ParameterHandler/CheckDataConsistency.cuh"

#include "Definitions/Matrices.cuh"

namespace uammd{
namespace structured{

    template<class PairType>
    class PairParameterHandler{

        public:

            using PairParameters       = typename PairType::PairParameters;
            using InputPairParameters  = typename PairType::InputPairParameters;

        private:

            std::shared_ptr<GlobalData>    gd;
            std::shared_ptr<ParticleGroup> pg;
            std::shared_ptr<ParticleData>  pd;

            std::shared_ptr<Types::TypesHandler> types;

            std::map<std::tuple<std::string,std::string,int>,PairParameters> pairParameters;

            int nPairTypes;

            bool isBatched;
            int nBatches;

            thrust::host_vector<PairParameters>   pairParameters_cpu;
            thrust::device_vector<PairParameters> pairParameters_gpu;

            void add(int ti, int tj, int batchId, InputPairParameters p){

                auto newp = PairType::processPairParameters(p);

                UAMMD_SET_3D_ROW_MAJOR(pairParameters_cpu,nBatches,nPairTypes,nPairTypes,batchId,ti,tj,newp);
                UAMMD_SET_3D_ROW_MAJOR(pairParameters_cpu,nBatches,nPairTypes,nPairTypes,batchId,tj,ti,newp);

                //UAMMD_SET_3D_COL_MAJOR(pairParameters_cpu,nBatches,nPairTypes,nPairTypes,batchId,ti,tj,newp);
                //UAMMD_SET_3D_COL_MAJOR(pairParameters_cpu,nBatches,nPairTypes,nPairTypes,batchId,tj,ti,newp);

                ////////////////////////////////////////

                std::tuple<std::string,std::string,int> key = std::make_tuple(p.name_i, p.name_j, batchId);
                //Check if the pair parameters have already been added
                if(pairParameters.find(key) == pairParameters.end()){
                    pairParameters[key] = newp;

                    //Check if the inverse pair parameters have already been added
                    std::tuple<std::string,std::string,int> key2 = std::make_tuple(p.name_j, p.name_i, batchId);
                    if(pairParameters.find(key2) != pairParameters.end() and p.name_i != p.name_j){
                        System::log<System::DEBUG>("[PairParameterHandler] Added pair paraemters (%s,%s) but inverse pair parameters (%s,%s) already exist. The inverse pair parameters will be overwritten.",
                                                    p.name_i.c_str(),p.name_j.c_str(),p.name_j.c_str(),p.name_i.c_str());
                    }
                }
                else{
                    System::log<System::CRITICAL>("[PairParameterHandler] Trying to add pair parameters for the same pair (%s,%s) twice!",
                                                  p.name_i.c_str(),p.name_j.c_str());
                }
            }


        public:

            PairParameterHandler(std::shared_ptr<GlobalData>    gd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 DataEntry& data):gd(gd), pg(pg),
                                                  pd(pg->getParticleData()),
                                                  types(gd->getTypes()){

                auto pairsData = data.getDataMap();

                ///////////////////////////////////////////////////////////////////

                //Check if "batchId" is in labels
                std::vector<std::string> labels = data.getLabels();
                if(std::find(labels.begin(), labels.end(), "batchId") == labels.end()){
                    System::log<System::MESSAGE>("[PairParameterHandler] Parameters are not batched");
                    nPairTypes = types->getNumberOfTypes();

                    isBatched = false;
                    nBatches  = 1;

                } else {
                    System::log<System::MESSAGE>("[PairParameterHandler] Parameters are batched");

                    if(!Batching::checkDataConsistency(data)){
                        System::log<System::CRITICAL>("[PairParameterHandler] Parameters are not consistent");
                    }

                    if(!Batching::checkParticleGroupDataConsistency(pg, data)){
                        System::log<System::CRITICAL>("[PairParameterHandler] Parameters are not consistent with particle group");
                    }

                    nPairTypes = types->getNumberOfTypes();

                    isBatched = true;
                    std::vector<int> batchIds = data.getData<int>("batchId");
                    //Count number of different batch ids. Convert to set to remove duplicates
                    nBatches = std::set<int>(batchIds.begin(), batchIds.end()).size();
                }

                System::log<System::MESSAGE>("[PairParameterHandler] Number of pairs types: %d", nPairTypes);
                System::log<System::MESSAGE>("[PairParameterHandler] Number of batches: %d", nBatches);

                pairParameters_cpu.resize(nBatches*nPairTypes*nPairTypes);

                ///////////////////////////////////////////////////////////////////

                InputPairParameters pairParamBuffer;
                for(auto pairInfo : pairsData){

                    pairParamBuffer = PairType::readPairParameters(pairInfo);

                    std::string name_i = pairParamBuffer.name_i;
                    std::string name_j = pairParamBuffer.name_j;

                    int ti = types->getTypeId(name_i);
                    int tj = types->getTypeId(name_j);

                    if(isBatched){
                        int batchId = pairInfo["batchId"];

                        this->add(ti,tj,batchId,pairParamBuffer);

                        ExtendedSystem::log<ExtendedSystem::DEBUG>("[PairParameterHandler] Pair %s,%s (type id: %i,%i, batch id: %i) added.",
                                                                    name_i.c_str(),name_j.c_str(),
                                                                    ti,tj,batchId);
                    } else {
                        this->add(ti,tj,0,pairParamBuffer);

                        ExtendedSystem::log<ExtendedSystem::DEBUG>("[PairParameterHandler] Pair %s,%s (type id: %i,%i) added.",
                                                                    name_i.c_str(),name_j.c_str(),
                                                                    ti,tj);
                    }
                }

                pairParameters_gpu = pairParameters_cpu;

            }

            ~PairParameterHandler(){}

            int getNumTypes()  { return nPairTypes; }
            int getNumBatches(){ return nBatches; }

            std::map<std::tuple<std::string,std::string,int>,PairParameters> getPairParameters(){
                return pairParameters;
            }

            struct PairIterator{

                real4* pos;
                int*   batchId;

                PairParameters* pairParam;

                int nPairTypes;
                int nBatches;

                PairIterator(){}
                PairIterator(real4* pos,
                             int*   batchId,
                             PairParameters* pairParam,
                             int nPairTypes,
                             int nBatches):pos(pos),batchId(batchId),
                                           pairParam(pairParam),
                                           nPairTypes(nPairTypes),
                                           nBatches(nBatches){}

                inline __device__ PairParameters operator()(const int& index_i,const int& index_j) const {

                    int ti = int(0);
                    int tj = int(0);
                    if(nPairTypes > 1){
                        ti = int(pos[index_i].w);
                        tj = int(pos[index_j].w);
                    }

                    int s = int(0);
                    if(nBatches > 1){
                        s = batchId[index_i]; //We assume that batchId is the same for both particles
                    }

                    cub::CacheModifiedInputIterator<cub::LOAD_CA, PairParameters> itr(pairParam);

                    return UAMMD_GET_3D_ROW_MAJOR(itr,nBatches,nPairTypes,nPairTypes,s,ti,tj);
                    //return UAMMD_GET_3D_COL_MAJOR(itr,nBatches,nPairTypes,nPairTypes,s,ti,tj);
                }
            };

            PairIterator getPairIterator(){
                auto pos   = pd->getPos(access::location::gpu, access::mode::read);
                auto batchId = pd->getBatchId(access::location::gpu, access::mode::read);

                auto pp = thrust::raw_pointer_cast(pairParameters_gpu.data());

                return PairIterator(pos.raw(),batchId.raw(),pp,nPairTypes,nBatches);
            }
    };
}}
