#ifndef __VERLET_CONDITIONAL_LIST_SET_BASE__
#define __VERLET_CONDITIONAL_LIST_SET_BASE__

#include"utils/Grid.cuh"
#include<thrust/device_vector.h>

namespace uammd{
namespace structured{

    //Interface
    class VerletConditionalListSetBase {

        protected:

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;
            std::shared_ptr<ParticleGroup>         pg;

        public:

            using StrideIterator = cub::CountingInputIterator<int>;

            struct NeighbourListData{
                int   N;
                const int* neighbourList;
                const int* numberNeighbours;
                StrideIterator neighbourStart = StrideIterator(0);
            };

            VerletConditionalListSetBase(std::shared_ptr<GlobalData>    gd,
                                         std::shared_ptr<ParticleGroup> pg):gd(gd),pd(getExtendedParticleData(pg->getParticleData())),pg(pg){

                System::log<System::MESSAGE>("[VerletConditionalListSetBase] Constructing VerletConditionalListSetBase");
            }

            virtual std::string getName() = 0;

            virtual void update(cudaStream_t st) = 0;
            virtual NeighbourListData getNeighbourList(std::string conditionName) = 0;

            virtual void setCutOff(real cutOff) = 0;
            virtual real getCutOff() = 0;
    };

}}


#endif
