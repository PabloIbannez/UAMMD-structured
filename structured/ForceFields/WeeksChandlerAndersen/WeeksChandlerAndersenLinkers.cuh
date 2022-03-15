#ifndef __WCA_LINKERS__
#define __WCA_LINKERS__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace WeeksChandlerAndersen{
    
    template<class Base_  >
    class WeeksChandlerAndersenLinkers : public Disordered::LinkersBonded<WeeksChandlerAndersen<Base_>>{

            using Base = Disordered::LinkersBonded<WeeksChandlerAndersen<Base_>>;

        protected:
            
            using NeighbourList = FixedNeighbourList<typename Base::Topology>;

            using InteractorWCALinkersType = Interactor::PairInteractor<typename Base::WCAType,NeighbourList>;

            std::shared_ptr<NeighbourList> fnl; 
            
            std::shared_ptr<InteractorWCALinkersType> wcalinkers;
        
        public:
        
            WeeksChandlerAndersenLinkers(std::shared_ptr<System>        sys,
                                         std::shared_ptr<ParticleData>  pd,
                                         std::shared_ptr<ParticleGroup> pg,
                                         InputFile&                     in):Base(sys,pd,pg,in){
                
                //Fixed Neighbour list
                
                typename NeighbourList::Parameters NeighbourListParam;
                
                NeighbourListParam.topologyLabel = "LINKER_UNBOUND";

                fnl = std::make_shared<NeighbourList>(this->sys,this->pd,this->pg,
                                                      this->top,
                                                      NeighbourListParam);
                
                //Add interactor
                typename InteractorWCALinkersType::Parameters interactorWCALinkersParameters;

                interactorWCALinkersParameters.name = "WeeksChandlerAndersenLinkers";
                interactorWCALinkersParameters.pot  = this->potWCA;
                interactorWCALinkersParameters.nl   = fnl;
                interactorWCALinkersParameters.conditionInteractionName = "";// Ignored in fixed neig list

                wcalinkers = std::make_shared<InteractorWCALinkersType>(this->sys,this->pd,this->pg,
                                                                      interactorWCALinkersParameters);


            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                this->wcalinkers->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                this->wcalinkers->updateBox(box);
            }
    };


}}}}


#endif
