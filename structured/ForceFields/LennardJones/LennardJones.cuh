#ifndef __LJ__
#define __LJ__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace LennardJones{

    template<class Base_>
    class LennardJones : public Base_{

            using Base = Base_;
        
        protected:
            
            using LJType  = Potentials::UnBound::LennardJones;
            
            using InteractorLJType = Interactor::PairInteractor<LJType,typename Base::NeighbourList>;
            
            std::shared_ptr<typename LJType::ParameterPairsHandler> ljparam;

            std::shared_ptr<LJType>        potLJ;
            std::shared_ptr<InteractorLJType> lj;
        
        
        protected:
            
            real cutOffDstLJ;

        public:
        
            LennardJones(std::shared_ptr<System>        sys,
                         std::shared_ptr<ParticleData>  pd,
                         std::shared_ptr<ParticleGroup> pg,
                         InputFile&                     in):Base(sys,pd,pg,in),
                                                            cutOffDstLJ(std::stof(in.getOption("cutOffDstLJ",InputFile::Required).str())){
                
                ljparam  = this->top->template readPairs<typename LJType::InteractionParameters>("LJ");
                
                if(cutOffDstLJ >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[LennardJones] cutOffDstLJ (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstLJ,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[LennarJones] "
                                                         "Parameter cutOffDstLJ added: %f",
                                                          cutOffDstLJ);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDstLJ));
                
                //Add lj
                typename LJType::Parameters ljPotentialParam;
                
                ljPotentialParam.paramPairsHandler = ljparam;
                ljPotentialParam.cutOffDstLJ = cutOffDstLJ;
                
                std::shared_ptr<LJType> potLJ = std::make_shared<LJType>(this->pd,ljPotentialParam);

                typename InteractorLJType::Parameters interactorLJParameters;

                interactorLJParameters.name = "LJ";
                interactorLJParameters.pot  = potLJ;
                interactorLJParameters.nl   = this->nl;
                interactorLJParameters.conditionInteractionName = "inter";

                lj = std::make_shared<InteractorLJType>(this->sys,this->pd,this->pg,
                                                        interactorLJParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                lj->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                lj->updateBox(box);
            }
    
    };

}}}}


#endif
