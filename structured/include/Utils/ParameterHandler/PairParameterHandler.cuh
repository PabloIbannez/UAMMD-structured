#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Utils/ParameterHandler/CheckDataConsistency.cuh"

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

            void add(int ti, int tj, int batchId, InputPairParameters p);


        public:

            PairParameterHandler(std::shared_ptr<GlobalData>    gd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 DataEntry& data);

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

            PairIterator getPairIterator();
    };
}}
