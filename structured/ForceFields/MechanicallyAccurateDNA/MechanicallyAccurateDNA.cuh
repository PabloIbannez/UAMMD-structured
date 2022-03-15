#ifndef __MADNA__
#define __MADNA__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace MechanicallyAccurateDNA{

    template<class Base_ >
    class MechanicallyAccurateDNA : public Base_{

            using Base = Base_;
        
        protected:

            const real epsilon = 1.0;
            
            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::HarmonicAngular;
            using DihedralType = Potentials::Bond4::Dihedral;

            using DHType    = Potentials::UnBound::DebyeHuckel<typename Base::Units>;
            using WCAType  = Potentials::UnBound::WCA;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            using InteractorDihedralType   = Interactor::BondedInteractor<DihedralType,
                                                                          Interactor::BondedInteractor_ns::BondProcessor<DihedralType>,
                                                                          Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,DihedralType>>;
            
            using InteractorDHType   = Interactor::PairInteractor<DHType,typename Base::NeighbourList>;
            using InteractorWCAType = Interactor::PairInteractor<WCAType,typename Base::NeighbourList>;
            

            std::shared_ptr<InteractorBondType>     bonds;
            std::shared_ptr<InteractorAngleType>    angles;
            std::shared_ptr<InteractorDihedralType> dihedrals;
            
            std::shared_ptr<InteractorDHType>  dh;
            std::shared_ptr<InteractorWCAType> wca;
        
        
        protected:
            
            real cutOffDstDH;
            real cutOffDstWCA;

            real dielectricConstant;
            real debyeLength;

        public:
        
            MechanicallyAccurateDNA(std::shared_ptr<System>        sys,
                                    std::shared_ptr<ParticleData>  pd,
                                    std::shared_ptr<ParticleGroup> pg,
                                    InputFile&                     in):Base(sys,pd,pg,in),
                                                                       cutOffDstDH(std::stof(in.getOption("cutOffDstDH",InputFile::Required).str())),
                                                                       cutOffDstWCA(std::stof(in.getOption("cutOffDstWCA",InputFile::Required).str())),
                                                                       dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                                       debyeLength(std::stof(in.getOption("debyeLength",InputFile::Required).str())){
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[MechanicallyAccurateDNA] Mechanically Accurate DNA force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }

                if(cutOffDstDH >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[MechanicallyAccurateDNA] cutOffDstDH (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstDH,this->nl->getCutOffVerlet());
                }
                
                if(cutOffDstWCA >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[MechanicallyAccurateDNA] cutOffDstWCA (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstWCA,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[MechanicallyAccurateDNA] "
                                                         "Parameter cutOffDstDH added: %f",
                                                          cutOffDstDH);
                this->sys->template log<System::MESSAGE>("[MechanicallyAccurateDNA] "
                                                         "Parameter cutOffDstWCA added: %f",
                                                          cutOffDstWCA);
                
                this->sys->template log<System::MESSAGE>("[MechanicallyAccurateDNA] "
                                                         "Parameter dielectricConstant added: %f",
                                                          dielectricConstant);
                this->sys->template log<System::MESSAGE>("[MechanicallyAccurateDNA] "
                                                         "Parameter debyeLength added: %f",
                                                          debyeLength);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             std::max(cutOffDstDH,cutOffDstWCA)));
                
                //Add bonds
                BondType::Parameters bondParameters;

                std::shared_ptr<BondType> b_PBC = std::make_shared<BondType>(this->pd,
                                                                             bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, b_PBC,
                                                             interactorBondParameters);
                
                //Add angles
                AngleType::Parameters angleParameters;

                std::shared_ptr<AngleType> a_PBC = std::make_shared<AngleType>(this->pd,
                                                                               angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                               this->top, a_PBC,
                                                               interactorAngleParameters);
                
                //Add dihedrals
                DihedralType::Parameters dihedralParameters;

                std::shared_ptr<DihedralType> d_PBC = std::make_shared<DihedralType>(this->pd,
                                                                                     dihedralParameters);
                
                typename InteractorDihedralType::Parameters interactorDihedralParameters;
                
                interactorDihedralParameters.bondName = "DIHEDRALS";

                dihedrals = std::make_shared<InteractorDihedralType>(this->sys, this->pd, this->pg,
                                                                     this->top, d_PBC,
                                                                     interactorDihedralParameters);
                
                //Add dh
                typename DHType::Parameters dhPotentialParam;
                
                dhPotentialParam.dielectricConstant = dielectricConstant;
                dhPotentialParam.debyeLength        = debyeLength;
                
                dhPotentialParam.cutOff = cutOffDstDH;
                
                std::shared_ptr<DHType> potDH = std::make_shared<DHType>(this->pd,dhPotentialParam);

                typename InteractorDHType::Parameters interactorDHParameters;

                interactorDHParameters.name = "DH";
                interactorDHParameters.pot  = potDH;
                interactorDHParameters.nl   = this->nl;
                interactorDHParameters.conditionInteractionName = "charged";

                dh = std::make_shared<InteractorDHType>(this->sys,this->pd,this->pg,
                                                        interactorDHParameters);
                
                //Add wca
                typename WCAType::Parameters wcaPotentialParam;
                
                wcaPotentialParam.epsilon = epsilon;
                
                wcaPotentialParam.cutOff = cutOffDstWCA;
                
                std::shared_ptr<WCAType> potWCA = std::make_shared<WCAType>(this->pd,wcaPotentialParam);

                typename InteractorWCAType::Parameters interactorWCAParameters;

                interactorWCAParameters.name = "WCA";
                interactorWCAParameters.pot  = potWCA;
                interactorWCAParameters.nl   = this->nl;
                interactorWCAParameters.conditionInteractionName = "nonExcluded";

                wca = std::make_shared<InteractorWCAType>(this->sys,this->pd,this->pg,
                                                          interactorWCAParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
                angles->sum(comp,st);
                dihedrals->sum(comp,st);
                dh->sum(comp,st);
                wca->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
                angles->updateBox(box);
                dihedrals->updateBox(box);
                dh->updateBox(box);
                wca->updateBox(box);
            }
    
    };

}}}}


#endif
