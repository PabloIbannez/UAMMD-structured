#ifndef __PAIR_INTERACTION__
#define __PAIR_INTERACTION__

#include"Interactor/NeighbourList/CellList.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace PairInteractor_ns{
    
    template<class PotentialTransverser, class NeighListData>
    __global__ void transverseNeighboursThreadPerParticle(int numberParticles,
                                                          const real4*  __restrict__ pos,
                                                          ParticleGroup::IndexIterator group2GlobalIndex,
                                                          NeighListData nlData,
                                                          PotentialTransverser potTransverser){
        int id = blockIdx.x*blockDim.x + threadIdx.x;

        if(id>=numberParticles){return;}

        const int nn   = nlData.numberNeighbours[id];
        
        if(nn>int(1)){
        
            int const * nl = nlData.neighbourList + nlData.neighbourStart[id];

            const int i_global = nl[0];
            const real4 pi     = pos[i_global];
            
            typename PotentialTransverser::resultType quantityLocal = potTransverser.zero();
            
            for(int n = int(1); n<nn; n += int(1)) {

                const int j_global = nl[n*numberParticles];
                const real4 pj     = pos[j_global];

                potTransverser.accumulate(quantityLocal,potTransverser.compute(i_global,
                                                                               j_global,
                                                                               pi,pj));
            }
            
            potTransverser.set(i_global,quantityLocal);
        }
    }
    
    }
    
    template<class Potential, class NeighbourList, int THREADS_PER_BLOCK = 256>
    class PairInteractor: public Interactor{
        
        std::shared_ptr<Potential>     pot;
        std::shared_ptr<NeighbourList> nl;

        Box box;

        std::string conditionInteractionName;

        public:

            struct Parameters{

                std::string name = "PairInteractor";
                
                std::shared_ptr<Potential>     pot;
                std::shared_ptr<NeighbourList> nl;
                
                std::string conditionInteractionName;

            };

            PairInteractor(std::shared_ptr<ParticleGroup> pg, 
                           Parameters par):Interactor(pg,par.name),
                                           pot(par.pot),
                                           nl(par.nl),
                                           conditionInteractionName(par.conditionInteractionName){}

            void sum(Computables comp,cudaStream_t st) override {
                
                nl->update(st);
                
                auto nlData   = nl->getNeighbourList(conditionInteractionName,st);

                auto pos      = pd->getPos(access::location::gpu, access::mode::read);
                
                int numberParticles    = nlData.N; //numberParticles in the neig list
                auto group2GlobalIndex = pg->getIndexIterator(access::location::gpu);

                int Nthreads=THREADS_PER_BLOCK;
                int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);
                
                if(comp.force == true){
                
                    PairInteractor_ns::transverseNeighboursThreadPerParticle
                    <typename Potential::forceTransverser,typename NeighbourList::NeighbourListData>
                    <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                 pos.raw(),
                                                 group2GlobalIndex,
                                                 nlData,
                                                 pot->getForceTransverser());
                    CudaCheckError();
                }

                if(comp.energy == true){

                    PairInteractor_ns::transverseNeighboursThreadPerParticle
                    <typename Potential::energyTransverser,typename NeighbourList::NeighbourListData>
                    <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                 pos.raw(),
                                                 group2GlobalIndex,
                                                 nlData,
                                                 pot->getEnergyTransverser());
                    CudaCheckError();
                }
                
                if(comp.virial == true){

                    PairInteractor_ns::transverseNeighboursThreadPerParticle
                    <typename Potential::virialTransverser,typename NeighbourList::NeighbourListData>
                    <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                 pos.raw(),
                                                 group2GlobalIndex,
                                                 nlData,
                                                 pot->getVirialTransverser());
                    CudaCheckError();
                }
            }
            
            void updateBox(Box newBox) override {
                box = newBox;
                pot->updateBox(box);
                nl->updateBox(box);
            }
    };

}}}

#endif
