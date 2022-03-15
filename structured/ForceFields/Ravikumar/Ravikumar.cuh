#ifndef __RAVIKUMAR__
#define __RAVIKUMAR__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Ravikumar{

    template<class Intra_  ,
             class RModel_ >
    class Ravikumar : public Intra_ {

            using Base = Intra_;
        
        protected:
            
            static constexpr real refTemperature = real(300.0);

            using RType    = Potentials::UnBound::RavikumarPotential<typename Base::Units>;
            using RModel   = RModel_;

            using InteractorRType   = Interactor::PairInteractor<RType,typename Base::NeighbourList>;
            
            using MJType   = Potentials::CommonParameters::StatisticalPotential::StatisticalPotential<typename Base::Topology>;

            std::shared_ptr<MJType> mj;
            
            std::shared_ptr<RType>         potR;
            std::shared_ptr<InteractorRType>  ravikumar;
            
            real d = real(3.8);
            real SASAThreshold = real(10.0);

        protected:
            
            real cutOffDstNP;
            real cutOffDstDH;

            real dielectricConstant;
            real debyeLenght;

        public:
        
            Ravikumar(std::shared_ptr<System>        sys,
                      std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<ParticleGroup> pg,
                      InputFile&                     in):Base(sys,pd,pg,in),
                                                         cutOffDstNP(std::stof(in.getOption("cutOffDstNP",InputFile::Required).str())),
                                                         cutOffDstDH(std::stof(in.getOption("cutOffDstDH",InputFile::Required).str())),
                                                         dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                         debyeLenght(std::stof(in.getOption("debyeLenght",InputFile::Required).str())){
                
                this->sys->template log<System::MESSAGE>("[Ravikumar] "
                                                         "Parameter refTemperature : %f",
                                                          refTemperature);
                
                if(cutOffDstNP >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[Ravikumar] cutOffDstNP (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstNP,this->nl->getCutOffVerlet());
                }

                if(cutOffDstDH >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[Ravikumar] cutOffDstDH (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstDH,this->nl->getCutOffVerlet());
                }

                this->nl->setCutOff(std::max(this->nl->getCutOff(),std::max(cutOffDstNP,cutOffDstNP)));
                
                this->sys->template log<System::MESSAGE>("[Ravikumar] "
                                                         "Parameter cutOffDstNP added: %f",
                                                          cutOffDstNP);
                this->sys->template log<System::MESSAGE>("[Ravikumar] "
                                                         "Parameter cutOffDstDH added: %f",
                                                          cutOffDstDH);
                
                this->sys->template log<System::MESSAGE>("[Ravikumar] "
                                                         "Parameter dielectricConstant added: %f",
                                                          dielectricConstant);
                this->sys->template log<System::MESSAGE>("[Ravikumar] "
                                                         "Parameter debyeLenght added: %f",
                                                          debyeLenght);
                
                mj = std::make_shared<MJType>(sys,pd,pg,this->top,refTemperature);
                mj->loadParameters("MJ",RModel::epsilon_0,RModel::lambda,true);

                //Load SASA

                this->top->template loadProperty(this->pd,"SASA",this->pd->getSASA(access::location::cpu, access::mode::write));

                //Set SASA threshold

                this->cond->setSASAThreshold(SASAThreshold);
                
                //Add Ravikumar
                typename RType::Parameters rPotentialParam;
                
                rPotentialParam.paramPairsHandler = mj->getParameters();
                
                rPotentialParam.dielectricConstant = dielectricConstant;
                rPotentialParam.debyeLenght        = debyeLenght;
                
                rPotentialParam.cutOffDstNP = cutOffDstNP;
                rPotentialParam.cutOffDstDH = cutOffDstDH;
                
                rPotentialParam.d     = d;
                rPotentialParam.gamma = RModel::gamma*real(2.0);
                
                potR = std::make_shared<RType>(this->pd,rPotentialParam);

                typename InteractorRType::Parameters interactorRParameters;

                interactorRParameters.name = "Ravikumar";
                interactorRParameters.pot  = potR;
                interactorRParameters.nl   = this->nl;
                interactorRParameters.conditionInteractionName = "inter";

                ravikumar = std::make_shared<InteractorRType>(this->sys,this->pd,this->pg,
                                                               interactorRParameters);

            }
            
            real getCutOff(std::string potName){
                if(potName == "Ravikumar"){
                    return std::max(cutOffDstNP,cutOffDstNP);
                }  

                this->sys->template log<System::CRITICAL>("[Ravikumar] Requested cutOff of the potential %s. "
                                                            "But %s is not present in the force field, or it has not got cutOff defined",
                                                            potName.c_str(),potName.c_str());
                return -1;

            }
            
            
            void sum(std::string potName,Computables comp,cudaStream_t st){
                if(potName == "Ravikumar"){
                    ravikumar->sum(comp,st);
                    return;
                }  

                this->sys->template log<System::CRITICAL>("[Ravikumar] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            potName.c_str(),potName.c_str());
            }

            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                ravikumar->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                ravikumar->updateBox(box);
            }
    
    };

}}}}


#endif
