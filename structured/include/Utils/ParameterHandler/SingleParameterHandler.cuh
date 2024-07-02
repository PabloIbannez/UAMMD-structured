#pragma once

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Utils/ParameterHandler/CheckDataConsistency.cuh"

namespace uammd{
namespace structured{

    template<class SingleType>
    class SingleParameterHandler{

        public:

            using SingleParameters      = typename SingleType::SingleParameters;
            using InputSingleParameters = typename SingleType::InputSingleParameters;

        private:

            std::shared_ptr<GlobalData>    gd;
            std::shared_ptr<ParticleGroup> pg;
            std::shared_ptr<ParticleData>  pd;

            std::shared_ptr<Types::TypesHandler> types;

            std::map<std::pair<std::string,int>, SingleParameters> singleParameters;

            int nSingleTypes;

            bool isBatched;
            int nBatches;

            thrust::host_vector<SingleParameters>   singleParameters_cpu;
            thrust::device_vector<SingleParameters> singleParameters_gpu;

            void add(int t, int batchId, InputSingleParameters p);


        public:

            SingleParameterHandler(std::shared_ptr<GlobalData>    gd,
                                   std::shared_ptr<ParticleGroup> pg,
                                   DataEntry& data);

            ~SingleParameterHandler(){}

            int getNumTypes()  { return nSingleTypes; }
            int getNumBatches(){ return nBatches; }

            std::map<std::pair<std::string,int>, SingleParameters> getSingleParameters(){
                return singleParameters;
            }

            struct SingleIterator{

                real4* pos;
                int*   batchId;

                SingleParameters* singleParam;

                int nSingleTypes;
                int nBatches;

                SingleIterator(){}
                SingleIterator(real4* pos,
                               int*   batchId,
                               SingleParameters * singleParam,
                               int nSingleTypes,
                               int nBatches):pos(pos),batchId(batchId),
                                             singleParam(singleParam),
                                             nSingleTypes(nSingleTypes),
                                             nBatches(nBatches){}

                inline __device__ SingleParameters operator()(const int& index) const{

                    int t = int(0);
                    if(nSingleTypes > 1){t = int(pos[index].w);}

                    int s = int(0);
                    if(nBatches > 1){s = batchId[index];}

                    cub::CacheModifiedInputIterator<cub::LOAD_CA, SingleParameters> itr(singleParam);
                    return UAMMD_GET_2D_ROW_MAJOR(itr, nBatches, nSingleTypes, s, t);
                    //return UAMMD_GET_2D_COL_MAJOR(itr, nBatches, nSingleTypes, s, t);
                }
            };

            SingleIterator getSingleIterator();

    };
}}

