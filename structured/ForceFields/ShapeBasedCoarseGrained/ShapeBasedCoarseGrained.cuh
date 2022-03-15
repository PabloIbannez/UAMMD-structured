#ifndef __SBCG__
#define __SBCG__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace ShapeBasedCoarseGrained{

    template<class Base_     >
    class ShapeBasedCoarseGrainedBase : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType    = Potentials::Bond2::Harmonic;
            using AngleType   = Potentials::Bond3::HarmonicAngular;
            using NativeType  = Potentials::Bond2::LennardJonesType3;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            using InteractorNativeType = Interactor::BondedInteractor<NativeType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<NativeType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,NativeType>>;
            
            std::shared_ptr<InteractorBondType>   bonds;
            //std::shared_ptr<InteractorAngleType>  angles;
            std::shared_ptr<InteractorNativeType> nativeContacts;

            bool bondsActive=true;
            bool nativeContactsActive=true;
        
        public:
        
            ShapeBasedCoarseGrainedBase(std::shared_ptr<System>        sys,
                                        std::shared_ptr<ParticleData>  pd,
                                        std::shared_ptr<ParticleGroup> pg,
                                        InputFile&                     in):Base(sys,pd,pg,in){
                                                                               
                //Add bonds
                BondType::Parameters bondParameters;
                std::shared_ptr<BondType> bp = std::make_shared<BondType>(this->pd,
                                                                          bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, bp,
                                                             interactorBondParameters);
                
                /*
                //Add angles
                AngleType::Parameters anglesParameters;
                std::shared_ptr<AngleType> ap = std::make_shared<AngleType>(this->pd,
                                                                            anglesParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                               this->top, ap,
                                                               interactorAngleParameters);*/

                if(this->top->isBlockPresent("PAIRS")){
                    //Add Native Contacts
                    NativeType::Parameters ncParameters;
                    std::shared_ptr<NativeType> NCp = std::make_shared<NativeType>(this->pd,
                                                                                   ncParameters);
                    
                    typename InteractorNativeType::Parameters interactorNativeContactParameters;
                    
                    interactorNativeContactParameters.bondName = "PAIRS";
                    
                    nativeContacts = std::make_shared<InteractorNativeType>(this->sys, this->pd, this->pg, 
                                                                            this->top, NCp,
                                                                            interactorNativeContactParameters);
                } else {
                    nativeContactsActive = false;
                    this->sys->template log<System::WARNING>("[ShapeBasedCoarseGrained] No PAIRS label");    
                }
                
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
                //angles->sum(comp,st);
                if(nativeContactsActive){nativeContacts->sum(comp,st);};
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
                //angles->updateBox(box);
                if(nativeContactsActive){nativeContacts->updateBox(box);};
            }
    
    };
    
    template<class Base_>
    class ShapeBasedCoarseGrained : public ShapeBasedCoarseGrainedBase<Base_>{
        
        protected:
            
            using Base = ShapeBasedCoarseGrainedBase<Base_>;
        
        protected:

            using StericType   = Potentials::UnBound::WCA;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,typename Base::NeighbourList>;
            
            std::shared_ptr<InteractorStericType> steric;
        
        protected:
            
            real cutOffDstSteric;
        
        public:
        
            ShapeBasedCoarseGrained(std::shared_ptr<System>        sys,
                                 std::shared_ptr<ParticleData>  pd,
                                 std::shared_ptr<ParticleGroup> pg,
                                 InputFile&                     in):Base(sys,pd,pg,in),
                                                                    cutOffDstSteric(std::stof(in.getOption("cutOffDstSteric",InputFile::Required).str())){
                
                if(cutOffDstSteric >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[ShapeBasedCoarseGrained] cutOffDstSteric (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstSteric,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[ShapeBasedCoarseGrained] "
                                                         "Parameter cutOffDstSteric added: %f",
                                                          cutOffDstSteric);    
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDstSteric));
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.cutOff = cutOffDstSteric;
                stericPotentialParam.epsilon = real(1.0);
                
                std::shared_ptr<StericType> pot = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "ShapeBasedCoarseGrainedSteric";
                interactorStericParameters.pot  = pot;
                interactorStericParameters.nl   = this->nl;
                interactorStericParameters.conditionInteractionName = "nonExcluded";

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
