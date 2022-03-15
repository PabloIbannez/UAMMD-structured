#ifndef __KB__
#define __KB__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace KaranicolasBrooks{

    template<class Base_     >
    class KaranicolasBrooksBase : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::HarmonicAngular;
            using DihedralType = Potentials::Bond4::Dihedral4;
            using NativeType   = Potentials::Bond2::LennardJonesKaranicolasBrooks;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;

            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;

            using InteractorDihedralType   = Interactor::BondedInteractor<DihedralType,
                                                                          Interactor::BondedInteractor_ns::BondProcessor<DihedralType>,
                                                                          Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,DihedralType>>;

            using InteractorNativeType = Interactor::BondedInteractor<NativeType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<NativeType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,NativeType>>;
            
            std::shared_ptr<InteractorBondType>     pairBonds;
            std::shared_ptr<InteractorAngleType>    angleBonds;
            std::shared_ptr<InteractorDihedralType> dihedralBonds;
            std::shared_ptr<InteractorNativeType>   nativeContacts;
        
        
        public:
        
            KaranicolasBrooksBase(std::shared_ptr<System>        sys,
                                  std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<ParticleGroup> pg,
                                  InputFile&                     in):Base(sys,pd,pg,in){
                
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[KaranicolasBrooks] Karanicolas Brooks force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }
                
                //Add pair bonds
                BondType::Parameters bondParameters;
                
                std::shared_ptr<BondType> bp_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "KB_BONDS";

                pairBonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                                 this->top, bp_PBC,
                                                                 interactorBondParameters);
                
                //Add angle bonds
                AngleType::Parameters angleParameters;
                
                std::shared_ptr<AngleType> ba_PBC = std::make_shared<AngleType>(this->pd,
                                                                                angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "KB_ANGLES";

                angleBonds = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                                   this->top, ba_PBC,
                                                                   interactorAngleParameters);
                
                //Add dihedral bonds
                DihedralType::Parameters dihedralParameters;
                
                std::shared_ptr<DihedralType> bd_PBC = std::make_shared<DihedralType>(this->pd,
                                                                                      dihedralParameters);
                
                typename InteractorDihedralType::Parameters interactorDihedralParameters;
                
                interactorDihedralParameters.bondName = "KB_DIHEDRALS";

                dihedralBonds = std::make_shared<InteractorDihedralType>(this->sys, this->pd, this->pg,
                                                                         this->top, bd_PBC,
                                                                         interactorDihedralParameters);
                
                //Add Native Contacts
                NativeType::Parameters ncParameters;
                
                std::shared_ptr<NativeType> NC_PBC = std::make_shared<NativeType>(this->pd,
                                                                                  ncParameters);
                
                typename InteractorNativeType::Parameters interactorNativeContactParameters;
                
                interactorNativeContactParameters.bondName = "KB_NATIVE_CONTACTS";
                
                nativeContacts = std::make_shared<InteractorNativeType>(this->sys, this->pd, this->pg, 
                                                                        this->top, NC_PBC,
                                                                        interactorNativeContactParameters);
                
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                pairBonds->sum(comp,st);
                angleBonds->sum(comp,st);
                dihedralBonds->sum(comp,st);
                nativeContacts->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                pairBonds->updateBox(box);
                angleBonds->updateBox(box);
                dihedralBonds->updateBox(box);
                nativeContacts->updateBox(box);
            }
    };
    
    template<class Base_>
    class KaranicolasBrooks : public KaranicolasBrooksBase<Base_>{
        
            using Base = KaranicolasBrooksBase<Base_>;
        
        protected:

            using StericType   = Potentials::UnBound::StericPBC<12>;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,typename Base::NeighbourList>;
            
            std::shared_ptr<InteractorStericType> steric;
        
        protected:
            
            real cutOffDstSteric;
        
        public:
        
            KaranicolasBrooks(std::shared_ptr<System>        sys,
                              std::shared_ptr<ParticleData>  pd,
                              std::shared_ptr<ParticleGroup> pg,
                              InputFile&                     in):Base(sys,pd,pg,in),
                                                                 cutOffDstSteric(std::stof(in.getOption("cutOffDstSteric",InputFile::Required).str())){
                
                if(cutOffDstSteric >= this->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[KaranicolasBrooks] cutOffDstSteric (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstSteric,this->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[KaranicolasBrooks] "
                                                         "Parameter cutOffDstSteric added: %f",
                                                          cutOffDstSteric);    
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDstSteric));

                //Load epsilon
                this->top->template loadProperty(this->pd,"EPSILON",this->pd->getEpsilon(access::location::cpu, access::mode::write));
                //Load innerRadius
                this->top->template loadProperty(this->pd,"INNER_RADIUS",this->pd->getInnerRadius(access::location::cpu, access::mode::write));
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.stericCutOff = cutOffDstSteric;
                
                stericPotentialParam.useInnerRadius = true;
                
                std::shared_ptr<StericType> pot = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "KaranicolasBrooksSteric";
                interactorStericParameters.pot  = pot;
                interactorStericParameters.nl   = this->nl;
                interactorStericParameters.conditionInteractionName = "intra";

                steric = std::make_shared<InteractorStericType>(this->sys,this->pd,this->pg,
                                                                interactorStericParameters);


                                                                    }
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                steric->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                steric->updateBox(box);
            }
    };

}}}}


#endif
