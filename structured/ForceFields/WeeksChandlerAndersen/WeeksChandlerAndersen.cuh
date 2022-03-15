#ifndef __WCA__
#define __WCA__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace WeeksChandlerAndersen{
    
    template<class Base_>
    class WeeksChandlerAndersen : public Base_{

            using Base = Base_;

        protected:
            
            using WCAType  = Potentials::UnBound::WCA;
            
            using InteractorWCAType = Interactor::PairInteractor<WCAType,typename Base::NeighbourList>;
            
            std::shared_ptr<WCAType>        potWCA;
            std::shared_ptr<InteractorWCAType> wca;
        
        
        protected:
            
            real epsilonWCA;
            real cutOffDstWCA;

        public:
        
            WeeksChandlerAndersen(std::shared_ptr<System>        sys,
                                  std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<ParticleGroup> pg,
                                  InputFile&                     in):Base(sys,pd,pg,in),
                                                                     cutOffDstWCA(std::stof(in.getOption("cutOffDstWCA",InputFile::Required).str())),
                                                                     epsilonWCA(std::stof(in.getOption("epsilonWCA",InputFile::Required).str())){
                
                if(cutOffDstWCA >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[WeeksChandlerAndersen] cutOffDstWCA (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstWCA,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[WeeksChandlerAndersen] "
                                                         "Parameter cutOffDstWCA added: %f",
                                                          cutOffDstWCA);
                
                this->sys->template log<System::MESSAGE>("[WeeksChandlerAndersen] "
                                                         "Parameter epsilonWCA added: %f",
                                                          epsilonWCA);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDstWCA));
                
                //Add wca
                typename WCAType::Parameters wcaPotentialParam;
                
                wcaPotentialParam.epsilon = epsilonWCA;
                
                wcaPotentialParam.cutOff = cutOffDstWCA;
                
                potWCA = std::make_shared<WCAType>(this->pd,wcaPotentialParam);

                typename InteractorWCAType::Parameters interactorWCAParameters;

                interactorWCAParameters.name = "WeeksChandlerAndersen";
                interactorWCAParameters.pot  = potWCA;
                interactorWCAParameters.nl   = this->nl;
                interactorWCAParameters.conditionInteractionName = "inter";

                wca = std::make_shared<InteractorWCAType>(this->sys,this->pd,this->pg,
                                                          interactorWCAParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                wca->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                wca->updateBox(box);
            }
    
    };
}}}}


#endif
