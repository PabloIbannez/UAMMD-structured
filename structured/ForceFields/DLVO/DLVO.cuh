#ifndef __DLVO__
#define __DLVO__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace DLVO{

    template<class Base_>
    class DLVO : public Base_{

            using Base = Base_;
        
        protected:

            //Non polar
            using NonPolarType           = Potentials::UnBound::LJ_WCA;
            using InteractorNonPolarType = Interactor::PairInteractor<NonPolarType,typename Base::NeighbourList>;
            
            std::shared_ptr<typename NonPolarType::ParameterPairsHandler> nonPolarParam;

            std::shared_ptr<NonPolarType>        potNonPolar;
            std::shared_ptr<InteractorNonPolarType> nonPolar;
                    
            bool nonPolarActive = false;

            //Polar
            using PolarType           = Potentials::UnBound::DebyeHuckelSpheres<typename Base::Units>;
            using InteractorPolarType = Interactor::PairInteractor<PolarType,typename Base::NeighbourList>;
    
            std::shared_ptr<PolarType>        potPolar;
            std::shared_ptr<InteractorPolarType> polar;
        
        protected:
            
            real cutOffDst;
            
            real dielectricConstant;
            real debyeLength;

        public:
        
            DLVO(std::shared_ptr<System>        sys,
                 std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<ParticleGroup> pg,
                 InputFile&                     in):Base(sys,pd,pg,in),
                                                    cutOffDst(std::stof(in.getOption("cutOffDst",InputFile::Required).str())),
                                                    dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                    debyeLength(std::stof(in.getOption("debyeLength",InputFile::Required).str())){
                
                
                if(cutOffDst >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[DLVO] cutOffDst (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDst,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[DLVO] "
                                                         "Parameter cutOffDst added: %f",
                                                          cutOffDst);
                
                this->sys->template log<System::MESSAGE>("[DLVO] "
                                                         "Parameter dielectricConstant added: %f",
                                                          dielectricConstant);
                this->sys->template log<System::MESSAGE>("[DLVO] "
                                                         "Parameter debyeLength added: %f",
                                                          debyeLength);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDst));
                
                //Add non polar
                if(this->top->isBlockPresent("NONPOLAR")){
                    nonPolarActive = true;
                
                    nonPolarParam  = this->top->template readPairs<typename NonPolarType::InteractionParameters>("NONPOLAR");
                    
                    typename NonPolarType::Parameters nonPolarPotentialParam;
                    
                    nonPolarPotentialParam.paramPairsHandler = nonPolarParam;
                    nonPolarPotentialParam.cutOffDst = cutOffDst;
                    
                    std::shared_ptr<NonPolarType> potNonPolar = std::make_shared<NonPolarType>(this->pd,nonPolarPotentialParam);

                    typename InteractorNonPolarType::Parameters interactorNonPolarParameters;

                    interactorNonPolarParameters.name = "NonPolar";
                    interactorNonPolarParameters.pot  = potNonPolar;
                    interactorNonPolarParameters.nl   = this->nl;
                    interactorNonPolarParameters.conditionInteractionName = "inter";

                    nonPolar = std::make_shared<InteractorNonPolarType>(this->sys,this->pd,this->pg,
                                                                        interactorNonPolarParameters);
                } else {
                    this->sys->template log<System::WARNING>("[DLVO] No NONPOLAR label");    
                }

                //Add polar
                typename PolarType::Parameters polarPotentialParam;

                polarPotentialParam.dielectricConstant = dielectricConstant;
                polarPotentialParam.debyeLenght        = debyeLength;
                polarPotentialParam.cutOff             = cutOffDst;
                    
                std::shared_ptr<PolarType> potPolar = std::make_shared<PolarType>(this->pd,polarPotentialParam);
    
                typename InteractorPolarType::Parameters interactorPolarParameters;

                interactorPolarParameters.name = "Polar";
                interactorPolarParameters.pot  = potPolar;
                interactorPolarParameters.nl   = this->nl;
                interactorPolarParameters.conditionInteractionName = "charged";

                polar = std::make_shared<InteractorPolarType>(this->sys,this->pd,this->pg,
                                                              interactorPolarParameters);
            }
            
            real getCutOff(std::string potName){
                if(potName == "DLVO"){
                    return cutOffDst;
                }  

                this->sys->template log<System::CRITICAL>("[DLVO] Requested cutOff of the potential %s. "
                                                            "But %s is not present in the force field, or it has not got cutOff defined",
                                                            potName.c_str(),potName.c_str());
                return -1;
            }
            
            void sum(std::string potName,Computables comp,cudaStream_t st){
                if(potName == "DLVO"){
                    if(nonPolarActive){nonPolar->sum(comp,st);}
                    polar->sum(comp,st);
                    return;
                }  

                this->sys->template log<System::CRITICAL>("[DLVO] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            potName.c_str(),potName.c_str());
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                if(nonPolarActive){nonPolar->sum(comp,st);}
                polar->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                if(nonPolarActive){nonPolar->updateBox(box);}
                polar->updateBox(box);
            }
    
    };

}}}}

#endif
