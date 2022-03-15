#ifndef __GENERIC__
#define __GENERIC__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Generic{

    template<class Base_     >
    class GenericBonded : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::HarmonicAngular;
            using DihedralType = Potentials::Bond4::Dihedral;
            using PairType     = Potentials::Bond2::LennardJonesType3;
            //using PairType     = Potentials::Bond2::LennardJonesGaussian;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;

            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;

            using InteractorDihedralType   = Interactor::BondedInteractor<DihedralType,
                                                                          Interactor::BondedInteractor_ns::BondProcessor<DihedralType>,
                                                                          Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,DihedralType>>;

            using InteractorPairType = Interactor::BondedInteractor<PairType,
                                                                    Interactor::BondedInteractor_ns::BondProcessor<PairType>,
                                                                    Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,PairType>>;
            
            std::shared_ptr<InteractorBondType>     pairBonds;
            std::shared_ptr<InteractorAngleType>    angleBonds;
            std::shared_ptr<InteractorDihedralType> dihedralBonds;
            std::shared_ptr<InteractorPairType>     pairs;
            
            bool bondsActive;
            bool anglesActive;
            bool dihedralsActive;
            bool pairsActive;
        
        public:
        
            GenericBonded(std::shared_ptr<System>        sys,
                          std::shared_ptr<ParticleData>  pd,
                          std::shared_ptr<ParticleGroup> pg,
                          InputFile&                     in):Base(sys,pd,pg,in){
                
                if(this->top->isBlockPresent("BONDS")){
                    //Add pair bonds
                    bondsActive = true;
                    BondType::Parameters bondParameters;
                    
                    std::shared_ptr<BondType> bp = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                    
                    typename InteractorBondType::Parameters interactorBondParameters;
                    
                    interactorBondParameters.bondName = "BONDS";

                    pairBonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                                     this->top, bp,
                                                                     interactorBondParameters);
                } else {
                    bondsActive = false;
                    this->sys->template log<System::WARNING>("[Generic] No BONDS label");    
                }
                
                if(this->top->isBlockPresent("ANGLES")){
                    //Add angle bonds
                    anglesActive = true;
                    AngleType::Parameters angleParameters;
                    
                    std::shared_ptr<AngleType> ba = std::make_shared<AngleType>(this->pd,
                                                                                angleParameters);
                    
                    typename InteractorAngleType::Parameters interactorAngleParameters;
                    
                    interactorAngleParameters.bondName = "ANGLES";

                    angleBonds = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                                       this->top, ba,
                                                                       interactorAngleParameters);
                } else {
                    anglesActive = false;
                    this->sys->template log<System::WARNING>("[Generic] No ANGLES label");    
                }
                
                if(this->top->isBlockPresent("DIHEDRALS")){
                    //Add dihedral bonds
                    dihedralsActive = true;
                    DihedralType::Parameters dihedralParameters;
                    
                    std::shared_ptr<DihedralType> bd = std::make_shared<DihedralType>(this->pd,
                                                                                      dihedralParameters);
                    
                    typename InteractorDihedralType::Parameters interactorDihedralParameters;
                    
                    interactorDihedralParameters.bondName = "DIHEDRALS";

                    dihedralBonds = std::make_shared<InteractorDihedralType>(this->sys, this->pd, this->pg,
                                                                             this->top, bd,
                                                                             interactorDihedralParameters);
                } else {
                    dihedralsActive = false;
                    this->sys->template log<System::WARNING>("[Generic] No DIHEDRALS label");    
                }
                
                if(this->top->isBlockPresent("PAIRS")){
                    //Add Pairs
                    pairsActive = true;
                    PairType::Parameters ncParameters;
                    
                    std::shared_ptr<PairType> pr = std::make_shared<PairType>(this->pd,
                                                                              ncParameters);
                    
                    typename InteractorPairType::Parameters interactorPairsParameters;
                    
                    interactorPairsParameters.bondName = "PAIRS";
                    
                    pairs = std::make_shared<InteractorPairType>(this->sys, this->pd, this->pg, 
                                                                 this->top, pr,
                                                                 interactorPairsParameters);
                } else {
                    pairsActive = false;
                    this->sys->template log<System::WARNING>("[Generic] No PAIRS label");    
                }
                
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                if(bondsActive)    {pairBonds->sum(comp,st);}
                if(anglesActive)   {angleBonds->sum(comp,st);}
                if(dihedralsActive){dihedralBonds->sum(comp,st);}
                if(pairsActive)    {pairs->sum(comp,st);}
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                if(bondsActive)    {pairBonds->updateBox(box);}
                if(anglesActive)   {angleBonds->updateBox(box);}
                if(dihedralsActive){dihedralBonds->updateBox(box);}
                if(pairsActive)    {pairs->updateBox(box);}
            }
    };
    
    template<class Base_>
    class GenericBondedWithSteric : public GenericBonded<Base_>{
        
        protected:
            
            using Base = GenericBonded<Base_>;
        
        protected:

            using StericType   = Potentials::UnBound::WCA;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,typename Base::NeighbourList>;
            
            std::shared_ptr<InteractorStericType> steric;
        
        protected:
            
            real cutOffDst;
        
        public:
        
            GenericBondedWithSteric(std::shared_ptr<System>        sys,
                                    std::shared_ptr<ParticleData>  pd,
                                    std::shared_ptr<ParticleGroup> pg,
                                    InputFile&                     in):Base(sys,pd,pg,in),
                                                                       cutOffDst(std::stof(in.getOption("cutOffDst",InputFile::Required).str())){
                
                if(cutOffDst >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[GenericBondedWithSteric] cutOffDst (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDst,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[GenericBondedWithSteric] "
                                                         "Parameter cutOffDst added: %f",
                                                          cutOffDst);    
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDst));
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.cutOff = cutOffDst;
                stericPotentialParam.epsilon = real(1.0);
                
                std::shared_ptr<StericType> pot = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "Steric";
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
