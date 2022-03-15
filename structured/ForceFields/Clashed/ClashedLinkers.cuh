#ifndef __CLASHED_LINKERS__
#define __CLASHED_LINKERS__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Clashed{

    template<class Intra_>
    class ClashedLinkers : public Disordered::LinkersBonded<Clashed<Intra_>>{

            using Base = Disordered::LinkersBonded<Clashed<Intra_>>;
        
        protected:
            
            using NeighbourList = FixedNeighbourList<typename Base::Topology>;

            using InteractorCLHDLinkersType = Interactor::PairInteractor<typename Base::CLHDType,NeighbourList>;

            std::shared_ptr<NeighbourList> nlf; 
            
            std::shared_ptr<InteractorCLHDLinkersType> clhdlinkers;
        
        public:
        
            ClashedLinkers(std::shared_ptr<System>        sys,
                           std::shared_ptr<ParticleData>  pd,
                           std::shared_ptr<ParticleGroup> pg,
                           InputFile&                     in):Base(sys,pd,pg,in){
                
                //Fixed Neighbour list
                
                typename NeighbourList::Parameters NeighbourListParam;
                
                NeighbourListParam.topologyLabel = "LINKER_UNBOUND";

                nlf = std::make_shared<NeighbourList>(this->sys,this->pd,this->pg,
                                                      this->top,
                                                      NeighbourListParam);
                
                //Add interactor
                typename InteractorCLHDLinkersType::Parameters interactorCLHDLinkersParameters;

                interactorCLHDLinkersParameters.name = "ClashedLinkers";
                interactorCLHDLinkersParameters.pot  = this->potCLHD;
                interactorCLHDLinkersParameters.nl   = nlf;
                interactorCLHDLinkersParameters.conditionInteractionName = "";// Ignored in fixed neig list

                clhdlinkers = std::make_shared<InteractorCLHDLinkersType>(this->sys,this->pd,this->pg,
                                                                          interactorCLHDLinkersParameters);


            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                this->clhdlinkers->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                this->clhdlinkers->updateBox(box);
            }
    };

}}}}


#endif
