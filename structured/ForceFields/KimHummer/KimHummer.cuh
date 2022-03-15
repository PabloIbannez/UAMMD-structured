#ifndef __KIM_HUMMER__
#define __KIM_HUMMER__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace KimHummer{

    template<class Base_  ,
             class KHModel_ >
    class KimHummer : public Base_ {

            using Base = Base_;
        
        protected:

            static constexpr real refTemperature = real(300.0);

            using KHType    = Potentials::UnBound::KimHummerPotential<typename Base::Units>;
            using KHModel   = KHModel_;

            using InteractorKHType   = Interactor::PairInteractor<KHType,typename Base::NeighbourList>;

            using MJType   = Potentials::CommonParameters::StatisticalPotential::StatisticalPotential<typename Base::Topology>;

            std::shared_ptr<MJType> mj;
            
            std::shared_ptr<KHType>         potKH;
            std::shared_ptr<InteractorKHType>  kh;

        protected:
            
            real cutOffDstNP;
            real cutOffDstDH;

            real dielectricConstant;
            real debyeLenght;

        public:
        
            KimHummer(std::shared_ptr<System>        sys,
                      std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<ParticleGroup> pg,
                      InputFile&                     in):Base(sys,pd,pg,in),
                                                         cutOffDstNP(std::stof(in.getOption("cutOffDstNP",InputFile::Required).str())),
                                                         cutOffDstDH(std::stof(in.getOption("cutOffDstDH",InputFile::Required).str())),
                                                         dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                         debyeLenght(std::stof(in.getOption("debyeLenght",InputFile::Required).str())){
                
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[KimHummer] Kim Hummer force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }
                
                this->sys->template log<System::MESSAGE>("[KimHummer] "
                                                         "Parameter refTemperature : %f",
                                                          refTemperature);
                
                if(cutOffDstNP >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[KimHummer] cutOffDstNP (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstNP,this->nl->getCutOffVerlet());
                }

                if(cutOffDstDH >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[KimHummer] cutOffDstDH (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstDH,this->nl->getCutOffVerlet());
                }

                this->nl->setCutOff(std::max(this->nl->getCutOff(),std::max(cutOffDstNP,cutOffDstNP)));
                
                this->sys->template log<System::MESSAGE>("[KimHummer] "
                                                         "Parameter cutOffDstNP added: %f",
                                                          cutOffDstNP);
                this->sys->template log<System::MESSAGE>("[KimHummer] "
                                                         "Parameter cutOffDstDH added: %f",
                                                          cutOffDstDH);
                
                this->sys->template log<System::MESSAGE>("[KimHummer] "
                                                         "Parameter dielectricConstant added: %f",
                                                          dielectricConstant);
                this->sys->template log<System::MESSAGE>("[KimHummer] "
                                                         "Parameter debyeLenght added: %f",
                                                          debyeLenght);
                
                mj = std::make_shared<MJType>(sys,pd,pg,this->top,refTemperature);
                mj->loadParameters("MJ",KHModel::epsilon_0,KHModel::lambda,true);

                //Load SASA

                if(KHModel::SASArequired){
                    this->top->template loadProperty(this->pd,"SASA",this->pd->getSASA(access::location::cpu, access::mode::write));
                    {
                    
                        auto typesparam = this->top->getTypes();

                        auto groupIndex = this->pg->getIndexIterator(access::location::cpu);
                        
                        auto pos  = this->pd->getPos(access::location::cpu, access::mode::readwrite);
                        auto SASA = this->pd->getSASA(access::location::cpu, access::mode::readwrite);
                        
                        fori(0,this->pg->getNumberParticles()){
                            int  index = groupIndex[i];
                            
                            real SASArc = typesparam->getTypeParameters(int(pos[index].w)).SASArc;
                            SASA[index]=SASA[index]/SASArc;
                            if(SASA[index] > real(1.0)){
                                SASA[index] = real(1.0);
                            }
                            SASA[index]=KHModel::SASAweight(SASA[index]/SASArc);
                        }
                            
                    }
                } else {
                        
                    auto groupIndex = this->pg->getIndexIterator(access::location::cpu);
                    
                    auto sasa = this->pd->getSASA(access::location::cpu, access::mode::readwrite);
                        
                    fori(0,this->pg->getNumberParticles()){
                        int  index = groupIndex[i];
                        sasa[index]=real(1.0);
                    }
                }
                
                //Add kh
                typename KHType::Parameters khPotentialParam;
                
                khPotentialParam.paramPairsHandler = mj->getParameters();
                
                khPotentialParam.dielectricConstant = dielectricConstant;
                khPotentialParam.debyeLenght        = debyeLenght;
                
                khPotentialParam.cutOffDstNP = cutOffDstNP;
                khPotentialParam.cutOffDstDH = cutOffDstDH;
                
                potKH = std::make_shared<KHType>(this->pd,khPotentialParam);

                typename InteractorKHType::Parameters interactorKHParameters;

                interactorKHParameters.name = "KH";
                interactorKHParameters.pot  = potKH;
                interactorKHParameters.nl   = this->nl;
                interactorKHParameters.conditionInteractionName = "inter";

                kh = std::make_shared<InteractorKHType>(this->sys,this->pd,this->pg,
                                                        interactorKHParameters);

            }
            
            real getCutOff(std::string potName){
                if(potName == "KimHummer"){
                    return std::max(cutOffDstNP,cutOffDstNP);
                }  

                this->sys->template log<System::CRITICAL>("[KimHummer] Requested cutOff of the potential %s. "
                                                            "But %s is not present in the force field, or it has not got cutOff defined",
                                                            potName.c_str(),potName.c_str());
                return -1;

            }
            
            
            void sum(std::string potName,Computables comp,cudaStream_t st){
                if(potName == "KimHummer"){
                    kh->sum(comp,st);
                    return;
                }  

                this->sys->template log<System::CRITICAL>("[KimHummer] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            potName.c_str(),potName.c_str());
            }

            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                kh->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                kh->updateBox(box);
            }
    
    };

}}}}


#endif
