#ifndef __SOP__
#define __SOP__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace SelfOrganizedPolymer{

    template<class Base_,class NativeType>
    class SelfOrganizedPolymerBase_ : public Base_{
        
        protected:

            using Base = Base_;
            
            static constexpr real K  = real(20.15);
            static constexpr real R0 = real(2.0);
            
            using BondType     = Potentials::Bond2::FeneConst_K_R0;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            using InteractorNativeType = Interactor::BondedInteractor<NativeType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<NativeType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,NativeType>>;
            
            std::shared_ptr<InteractorBondType>   covalentBonds;
            std::shared_ptr<InteractorNativeType> nativeContacts;
        
        
        public:
        
            SelfOrganizedPolymerBase_(std::shared_ptr<System>        sys,
                                      std::shared_ptr<ParticleData>  pd,
                                      std::shared_ptr<ParticleGroup> pg,
                                      InputFile&                     in):Base(sys,pd,pg,in){

                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[SelfOrganizedPolymer] Self organized polymer force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }
                
                //Add covalent bonds
                BondType::Parameters bondParameters;

                bondParameters.K   = K;
                bondParameters.R0  = R0;
                
                std::shared_ptr<BondType> b = std::make_shared<BondType>(this->pd,
                                                                         bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "SOP_BONDS";

                covalentBonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                                     this->top, b,
                                                                     interactorBondParameters);
                //Add Native Contacts
                typename NativeType::Parameters ncParameters;
                
                std::shared_ptr<NativeType> nc = std::make_shared<NativeType>(this->pd,
                                                                              ncParameters);
                
                typename InteractorNativeType::Parameters interactorNativeContactParameters;
                
                interactorNativeContactParameters.bondName = "SOP_NATIVE_CONTACT";
                
                nativeContacts = std::make_shared<InteractorNativeType>(this->sys, this->pd, this->pg, 
                                                                        this->top, nc,
                                                                        interactorNativeContactParameters);
                
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                covalentBonds->sum(comp,st);
                nativeContacts->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                covalentBonds->updateBox(box);
                nativeContacts->updateBox(box);
            }
    
    };

    template<class Base_,class NativeType_>
    class SelfOrganizedPolymer_ : public SelfOrganizedPolymerBase_<Base_,NativeType_>{
        
        protected:
            
            using Base = SelfOrganizedPolymerBase_<Base_,NativeType_>;

        protected:

            static constexpr real epsilon = real(1.0);
            static constexpr real sigma   = real(3.8);
            
            using StericType   = Potentials::UnBound::StericConstantPBC<6>;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,typename Base::NeighbourList>;
            
            std::shared_ptr<InteractorStericType> steric;
        
        protected:
            
            real cutOffDstSteric;
        
        public:
        
            SelfOrganizedPolymer_(std::shared_ptr<System>        sys,
                                  std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<ParticleGroup> pg,
                                  InputFile&                     in):Base(sys,pd,pg,in),
                                                                     cutOffDstSteric(std::stof(in.getOption("cutOffDstSteric",InputFile::Required).str())){
                
                if(cutOffDstSteric >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[SelfOrganizedPolymer] cutOffDstSteric (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                  cutOffDstSteric,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[SelfOrganizedPolymer] "
                                                         "Parameter cutOffDstSteric added: %f",
                                                          cutOffDstSteric);    

                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             cutOffDstSteric));
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.stericCutOff = cutOffDstSteric;
                
                stericPotentialParam.epsilon = epsilon;
                stericPotentialParam.sigma   = sigma;
                
                std::shared_ptr<StericType> pot = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "SelfOrganizedPolymerSteric";
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
    
    template<class Base_,class NativeType>
    class SelfOrganizedPolymerFixedSteric_: public SelfOrganizedPolymerBase_<Base_,NativeType>{
        
        protected:
            
            using Base = SelfOrganizedPolymerBase_<Base_,NativeType>;

            static constexpr real epsilon = real(1.0);
            static constexpr real sigma   = real(3.8);
            
            using FixedNeighbourList = FixedNeighbourList<typename Base::Topology>;
            
            using StericType    = Potentials::UnBound::StericConstantPBC<6>;
            
            using InteractorStericType = Interactor::PairInteractor<StericType,FixedNeighbourList>;
            
            std::shared_ptr<FixedNeighbourList> fnl; 
            
            std::shared_ptr<InteractorStericType> steric;
        
        protected:
            
            real cutOffDstSteric;
        
        public:
        
            SelfOrganizedPolymerFixedSteric_(std::shared_ptr<System>        sys,
                                             std::shared_ptr<ParticleData>  pd,
                                             std::shared_ptr<ParticleGroup> pg,
                                             InputFile&                     in):Base(sys,pd,pg,in),
                                                                                cutOffDstSteric(std::stof(in.getOption("cutOffDstSteric",InputFile::Required).str())){
                
                this->sys->template log<System::MESSAGE>("[SelfOrganizedPolymerFixedSteric] "
                                                         "Parameter cutOffDstSteric added: %f",
                                                          cutOffDstSteric);    
                //Fixed Neighbour list
                
                typename FixedNeighbourList::Parameters FixedNeighbourListParam;
                
                FixedNeighbourListParam.topologyLabel = "SOP_STERIC";

                fnl = std::make_shared<FixedNeighbourList>(this->sys,this->pd,this->pg,
                                                           this->top,
                                                           FixedNeighbourListParam);
                
                //Add steric potential
                StericType::Parameters stericPotentialParam;
                
                stericPotentialParam.stericCutOff = cutOffDstSteric;
                
                stericPotentialParam.epsilon = epsilon;
                stericPotentialParam.sigma   = sigma;
                
                std::shared_ptr<StericType> pot = std::make_shared<StericType>(this->pd,stericPotentialParam);

                typename InteractorStericType::Parameters interactorStericParameters;

                interactorStericParameters.name = "SelfOrganizedPolymerSteric";
                interactorStericParameters.pot  = pot;
                interactorStericParameters.nl   = fnl;
                interactorStericParameters.conditionInteractionName = ""; // Ignored in fixed neig list

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

    template<class Base>
    using SelfOrganizedPolymer = SelfOrganizedPolymer_<Base,Potentials::Bond2::LennardJonesType2>;
    
    template<class Base>
    using SelfOrganizedPolymerFixedSteric = SelfOrganizedPolymerFixedSteric_<Base,Potentials::Bond2::LennardJonesType2>;
    
    template<class Base>
    using SelfOrganizedPolymerGaussian = SelfOrganizedPolymer_<Base,Potentials::Bond2::LennardJonesGaussian>;
        

}}}}


#endif
