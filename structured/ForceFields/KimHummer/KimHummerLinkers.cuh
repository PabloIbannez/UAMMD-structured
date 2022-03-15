#ifndef __KIM_HUMMER_LINKERS__
#define __KIM_HUMMER_LINKERS__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace KimHummer{

    template<class Intra_  ,
             class KHModel_ >
    class KimHummerLinkers : public Disordered::LinkersBonded<KimHummer<Intra_,KHModel_>>{

            using Base = Disordered::LinkersBonded<KimHummer<Intra_,KHModel_>>;
        
        protected:
            
            using NeighbourList = FixedNeighbourList<typename Base::Topology>;

            using InteractorKHLinkersType = Interactor::PairInteractor<typename Base::KHType,NeighbourList>;

            std::shared_ptr<NeighbourList> nlf; 
            
            std::shared_ptr<InteractorKHLinkersType> khlinkers;
        
        public:
        
            KimHummerLinkers(std::shared_ptr<System>        sys,
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
                typename InteractorKHLinkersType::Parameters interactorKHLinkersParameters;

                interactorKHLinkersParameters.name = "KimHummerLinkers";
                interactorKHLinkersParameters.pot  = this->potKH;
                interactorKHLinkersParameters.nl   = nlf;
                interactorKHLinkersParameters.conditionInteractionName = "";// Ignored in fixed neig list

                khlinkers = std::make_shared<InteractorKHLinkersType>(this->sys,this->pd,this->pg,
                                                                      interactorKHLinkersParameters);


            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                this->khlinkers->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                this->khlinkers->updateBox(box);
            }
    
    };

}}}}


#endif
