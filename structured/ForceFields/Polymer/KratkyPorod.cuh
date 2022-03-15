#ifndef __KRATKY_POROD__
#define __KRATKY_POROD__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace KratkyPorodModel{

    template<class Base_ >
    class KratkyPorodModel : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::HarmonicConst_K_r0;
            using AngleType    = Potentials::Bond3::KratkyPorodConst_K;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            
            using ConstraintBaseType = Constraint::ShakeBond<Constraint::ShakeBasic_>;

            using ConstraintType     = Constraint::SHAKE<ConstraintBaseType,
                                                         Interactor::BondedInteractor_ns::BondProcessor<ConstraintBaseType>,
                                                         Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,ConstraintBaseType>>;

            
            std::shared_ptr<InteractorBondType>   bonds;
            std::shared_ptr<InteractorAngleType>  angles;

            bool constrained = false;
            std::shared_ptr<ConstraintType>  constraint;
        
        protected:

            real b;
            real B;
            real K;
        
        public:

            bool isConstrained(){
                return constrained;
            }

            std::shared_ptr<ConstraintType> getConstraint(){
                return constraint;
            };

        
            KratkyPorodModel(std::shared_ptr<System>        sys,
                             std::shared_ptr<ParticleData>  pd,
                             std::shared_ptr<ParticleGroup> pg,
                             InputFile&                     in):Base(sys,pd,pg,in),
                                                                b(std::stof(in.getOption("b",InputFile::Required).str())),
                                                                B(std::stof(in.getOption("B",InputFile::Required).str())),
                                                                K(std::stof(in.getOption("K",InputFile::Required).str())){
                
                this->sys->template log<System::MESSAGE>("[KratkyPorodModel] "
                                                         "Parameter b added: %f",
                                                          b);
                
                this->sys->template log<System::MESSAGE>("[KratkyPorodModel] "
                                                         "Parameter B added: %f",
                                                          B);
                
                this->sys->template log<System::MESSAGE>("[KratkyPorodModel] "
                                                         "Parameter K added: %f",
                                                          K);
                
                //Add bonds
                BondType::Parameters bondParameters;

                bondParameters.K  = K;
                bondParameters.r0 = b;
                
                std::shared_ptr<BondType> bH_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, bH_PBC,
                                                             interactorBondParameters);

                //Constrained
                
                ConstraintBaseType::Parameters constraintBaseParameters;
                
                std::shared_ptr<ConstraintBaseType> cnstr = std::make_shared<ConstraintBaseType>(this->pd,
                                                                                                 constraintBaseParameters);
                
                typename ConstraintType::Parameters constraintParameters;
                
                constraintParameters.bondName = "BONDS";

                constraint = std::make_shared<ConstraintType>(this->sys, this->pd, this->pg,
                                                              this->top, cnstr,
                                                              constraintParameters);




                
                //Add angles
                AngleType::Parameters angleParameters;

                angleParameters.K = B;
                
                std::shared_ptr<AngleType> aKP_PBC = std::make_shared<AngleType>(this->pd,
                                                                                 angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                               this->top, aKP_PBC,
                                                               interactorAngleParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
                angles->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
                angles->updateBox(box);
            }
    
    };
    
    template<class Base_ >
    class KratkyPorodModelFixed : public KratkyPorodModel<Base_>{
        
        protected:

            using KP_Base = KratkyPorodModel<Base_>;
            
            using FixedType     = Potentials::Bond1::HarmonicConstCommon_K_r0;
            
            using InteractorFixedType   = Interactor::BondedInteractor<FixedType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<FixedType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename KP_Base::Base::Topology,FixedType>>;
            std::shared_ptr<InteractorFixedType> fixed;
        
        protected:

            real Kfixed;
        
        public:
        
            KratkyPorodModelFixed(std::shared_ptr<System>        sys,
                                  std::shared_ptr<ParticleData>  pd,
                                  std::shared_ptr<ParticleGroup> pg,
                                  InputFile&                     in):KP_Base(sys,pd,pg,in),
                                                                     Kfixed(std::stof(in.getOption("Kfixed",InputFile::Required).str())){
                
                this->sys->template log<System::MESSAGE>("[KratkyPorodModelFixed] "
                                                         "Parameter Kfixed added: %f",
                                                          Kfixed);
                
                //Add fixed
                FixedType::Parameters fixedParameters;

                fixedParameters.K  = Kfixed;
                fixedParameters.r0 = real(0.0);
                
                std::shared_ptr<FixedType> fH_PBC = std::make_shared<FixedType>(this->pd,
                                                                                fixedParameters);
                
                typename InteractorFixedType::Parameters interactorFixedParameters;
                
                interactorFixedParameters.bondName = "FIXED";

                fixed = std::make_shared<InteractorFixedType>(this->sys, this->pd, this->pg,
                                                              this->top, fH_PBC,
                                                              interactorFixedParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                KP_Base::sum(comp,st);
                fixed->sum(comp,st);
            }
            
            void updateBox(Box box){
                KP_Base::updateBox(box);
                fixed->updateBox(box);
            }
    
    };
    
}}}}


#endif
