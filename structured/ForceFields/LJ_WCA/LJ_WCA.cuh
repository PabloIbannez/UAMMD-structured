#ifndef __LJ_WCA__
#define __LJ_WCA__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace LJ_WCA{

    template<class Base_>
    class LJ_WCA : public Base_{

            using Base = Base_;
        
        protected:
            
            using LJ_WCAType  = Potentials::UnBound::LJ_WCA;
            
            using InteractorLJ_WCAType = Interactor::PairInteractor<LJ_WCAType,typename Base::NeighbourList>;
            
            std::shared_ptr<typename LJ_WCAType::ParameterPairsHandler> lj_wca_param;

            std::shared_ptr<LJ_WCAType>        potLJ_WCA;
            std::shared_ptr<InteractorLJ_WCAType> lj_wca;
        
        
        protected:
            
            real cutOffDst;

        public:
        
            LJ_WCA(std::shared_ptr<System>        sys,
                         std::shared_ptr<ParticleData>  pd,
                         std::shared_ptr<ParticleGroup> pg,
                         InputFile&                     in):Base(sys,pd,pg,in),
                                                            cutOffDst(std::stof(in.getOption("cutOffDst",InputFile::Required).str())){
                
                lj_wca_param  = this->top->template readPairs<typename LJ_WCAType::InteractionParameters>("LJ_WCA");
                
                if(cutOffDst >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[LJ_WCA] cutOffDst (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDst,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[LennarJones] "
                                                         "Parameter cutOffDst added: %f",
                                                          cutOffDst);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDst));
                
                //Add lj wca
                typename LJ_WCAType::Parameters lj_wca_PotentialParam;
                
                lj_wca_PotentialParam.paramPairsHandler = lj_wca_param;
                lj_wca_PotentialParam.cutOffDst = cutOffDst;
                
                std::shared_ptr<LJ_WCAType> potLJ_WCA = std::make_shared<LJ_WCAType>(this->pd,lj_wca_PotentialParam);

                typename InteractorLJ_WCAType::Parameters interactorLJ_WCAParameters;

                interactorLJ_WCAParameters.name = "LJ_WCA";
                interactorLJ_WCAParameters.pot  = potLJ_WCA;
                interactorLJ_WCAParameters.nl   = this->nl;
                interactorLJ_WCAParameters.conditionInteractionName = "nonExcluded";

                lj_wca = std::make_shared<InteractorLJ_WCAType>(this->sys,this->pd,this->pg,
                                                                interactorLJ_WCAParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                lj_wca->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                lj_wca->updateBox(box);
            }
    
    };

}}}}


#endif
