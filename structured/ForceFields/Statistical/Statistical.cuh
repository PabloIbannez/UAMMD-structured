#ifndef __STATISTICAL_FF__
#define __STATISTICAL_FF__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Statistical{

    template<class Base_>
    class Statistical : public Base_ {

            using Base = Base_;
        
        protected:

            static constexpr real refTemperature = real(300.0);
            
            static constexpr real epsilon = real(1.0);
            static constexpr real sigma   = real(4.0);

            using StatParamType   = Potentials::CommonParameters::StatisticalPotential::StatisticalPotential<typename Base::Topology>;
            
            using StericType      = Potentials::UnBound::StericConstantPBC<12>;
            using StatisticalType = Potentials::UnBound::StatisticalPotential;
            using DHType          = Potentials::UnBound::DebyeHuckel<typename Base::Units>;
            //using DHType          = Potentials::UnBound::DebyeHuckelDistanceDependentDielectric<typename Base::Units>;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,typename Base::NeighbourList>;
            using InteractorDHType          = Interactor::PairInteractor<DHType,typename Base::NeighbourList>;
            using InteractorStatisticalType = Interactor::PairInteractor<StatisticalType,typename Base::NeighbourList>;
            
            std::shared_ptr<StatParamType> statParam;
            
            std::shared_ptr<InteractorStericType>     steric;
            std::shared_ptr<InteractorDHType>             dh;
            std::shared_ptr<InteractorStatisticalType>  stat;

        protected:
            
            real cutOffDstSteric;
            real cutOffDstDH;
            real cutOffDstStat;

            real dielectricConstant;
            real debyeLenght;

            real nContact;
            real r0Contact;

        public:
        
            Statistical(std::shared_ptr<System>        sys,
                        std::shared_ptr<ParticleData>  pd,
                        std::shared_ptr<ParticleGroup> pg,
                        InputFile&                     in):Base(sys,pd,pg,in),
                                                            cutOffDstSteric(std::stof(in.getOption("cutOffDstSteric",InputFile::Required).str())),
                                                           cutOffDstDH(std::stof(in.getOption("cutOffDstDH",InputFile::Required).str())),
                                                           cutOffDstStat(std::stof(in.getOption("cutOffDstStat",InputFile::Required).str())),
                                                           dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                           debyeLenght(std::stof(in.getOption("debyeLenght",InputFile::Required).str())),
                                                           nContact(std::stof(in.getOption("nContact",InputFile::Required).str())),
                                                           r0Contact(std::stof(in.getOption("r0Contact",InputFile::Required).str())){
                
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[Statistical] Statistical force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter refTemperature : %f",
                                                          refTemperature);
                
                if(cutOffDstSteric >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[Statistical] cutOffDstSteric (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstSteric,this->nl->getCutOffVerlet());
                }
                
                if(cutOffDstStat >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[Statistical] cutOffDstStat (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstStat,this->nl->getCutOffVerlet());
                }
                
                if(cutOffDstDH >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[Statistical] cutOffDstDH (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstDH,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter cutOffDstSteric added: %f",
                                                          cutOffDstSteric);
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter cutOffDstDH added: %f",
                                                          cutOffDstDH);
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter cutOffDstStat added: %f",
                                                          cutOffDstStat);
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter dielectricConstant added: %f",
                                                          dielectricConstant);
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter debyeLenght added: %f",
                                                          debyeLenght);
                
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter nContact added: %f",
                                                          nContact);
                this->sys->template log<System::MESSAGE>("[Statistical] "
                                                         "Parameter r0Contact added: %f",
                                                          r0Contact);

                this->nl->setCutOff(std::max(std::max(this->nl->getCutOff(),cutOffDstStat),cutOffDstDH));
                
                //

                auto& statLabelsStream = in.getOption("statLabels",InputFile::Required);

                std::vector<std::string> statLabels;
                std::string label;
                while(statLabelsStream >> label){
                    statLabels.push_back(label);
                }
                
                std::vector<real> epsilon_0;
                std::vector<real> lambda;

                auto& epsilon_0Stream = in.getOption("epsilon_0",InputFile::Required);
                real eps_0;
                while(epsilon_0Stream >> eps_0){epsilon_0.push_back(eps_0);}
                auto& lambdaStream = in.getOption("lambda",InputFile::Required);
                real lmd;
                while(lambdaStream >> lmd){lambda.push_back(lmd);}

                if(epsilon_0.size() != statLabels.size()){
                    this->sys->template log<System::CRITICAL>("[Statistical] The number of epsilon_0 parameters (%i),"
                                                              " does not match with the number of given labels (%i).",
                                                              epsilon_0.size(),statLabels.size());
                }
                
                if(lambda.size() != statLabels.size()){
                    this->sys->template log<System::CRITICAL>("[Statistical] The number of lambda parameters (%i),"
                                                              " does not match with the number of given labels (%i).",
                                                              lambda.size(),statLabels.size());
                }

                //for(uint i=0; i<statLabels.size(); i++){
                //    std::cout << statLabels[i] << " "
                //              << epsilon_0[i]  << " "
                //              << lambda[i]     << std::endl;
                //}

                //

                statParam = std::make_shared<StatParamType>(sys,pd,pg,this->top,refTemperature);
                for(uint i=0; i<statLabels.size(); i++){
                    statParam->loadParameters(statLabels[i],epsilon_0[i],lambda[i],false);
                }
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.stericCutOff = cutOffDstSteric;
                
                stericPotentialParam.epsilon = epsilon;
                stericPotentialParam.sigma   = sigma;
                
                std::shared_ptr<StericType> potSteric = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "Steric";
                interactorStericParameters.pot  = potSteric;
                interactorStericParameters.nl   = this->nl;
                interactorStericParameters.conditionInteractionName = "intra";

                steric = std::make_shared<InteractorStericType>(this->sys,this->pd,this->pg,
                                                                interactorStericParameters);

                //Add statistical
                typename StatisticalType::Parameters statPotentialParam;
                
                statPotentialParam.paramPairsHandler = statParam->getParameters();
                statPotentialParam.n                 = nContact;
                statPotentialParam.r0                = r0Contact;
                statPotentialParam.zeroEnergy        = real(0.1);
                statPotentialParam.cutOffDst         = cutOffDstStat;
                
                std::shared_ptr<StatisticalType> potStat = std::make_shared<StatisticalType>(this->pd,statPotentialParam);

                typename InteractorStatisticalType::Parameters interactorStatParameters;

                interactorStatParameters.name = "Statatistical";
                interactorStatParameters.pot  = potStat;
                interactorStatParameters.nl   = this->nl;
                interactorStatParameters.conditionInteractionName = "inter";

                stat = std::make_shared<InteractorStatisticalType>(this->sys,this->pd,this->pg,
                                                                   interactorStatParameters);
                
                //Add Debye-Huckel
                typename DHType::Parameters dhPotentialParam;
                
                dhPotentialParam.dielectricConstant = dielectricConstant;
                dhPotentialParam.debyeLenght        = debyeLenght;
                
                dhPotentialParam.cutOff = cutOffDstDH;
                
                std::shared_ptr<DHType> potDH = std::make_shared<DHType>(this->pd,dhPotentialParam);

                typename InteractorDHType::Parameters interactorDHParameters;

                interactorDHParameters.name = "DH";
                interactorDHParameters.pot  = potDH;
                interactorDHParameters.nl   = this->nl;
                interactorDHParameters.conditionInteractionName = "charged";

                dh = std::make_shared<InteractorDHType>(this->sys,this->pd,this->pg,
                                                        interactorDHParameters);

            }
            
            real getCutOff(std::string potName){
                if(potName == "Steric"){
                    return cutOffDstSteric;
                }  
                if(potName == "Statistical"){
                    return cutOffDstStat;
                }  
                if(potName == "DH"){
                    return cutOffDstDH;
                }  

                this->sys->template log<System::CRITICAL>("[Statistical] Requested cutOff of the potential %s. "
                                                            "But %s is not present in the force field, or it has not got cutOff defined",
                                                            potName.c_str(),potName.c_str());
                return -1;

            }
            
            
            void sum(std::string potName,Computables comp,cudaStream_t st){
                if(potName == "Steric"){
                    steric->sum(comp,st);
                    return;
                }  
                if(potName == "Statistical"){
                    stat->sum(comp,st);
                    return;
                }  
                if(potName == "DH"){
                    dh->sum(comp,st);
                    return;
                }  

                this->sys->template log<System::CRITICAL>("[Statatistical] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            potName.c_str(),potName.c_str());
            }

            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                steric->sum(comp,st);
                stat->sum(comp,st);
                dh->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                steric->updateBox(box);
                stat->updateBox(box);
                dh->updateBox(box);
            }
    
    };

}}}}


#endif
