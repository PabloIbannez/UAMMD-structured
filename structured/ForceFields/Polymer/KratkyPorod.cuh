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
            
            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::KratkyPorod;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            
            //using ConstraintBaseType = Constraint::ShakeBond<Constraint::ShakeBasic_>;

            //using ConstraintType     = Constraint::SHAKE<ConstraintBaseType,
            //                                             Interactor::BondedInteractor_ns::BondProcessor<ConstraintBaseType>,
            //                                             Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,ConstraintBaseType>>;

            
            std::shared_ptr<InteractorBondType>   bonds;
            std::shared_ptr<InteractorAngleType>  angles;

            //bool constrained = false;
            //std::shared_ptr<ConstraintType>  constraint;
        
        public:

            //bool isConstrained(){
            //    return constrained;
            //}

            //std::shared_ptr<ConstraintType> getConstraint(){
            //    return constraint;
            //};
        
            KratkyPorodModel(std::shared_ptr<System>        sys,
                             std::shared_ptr<ParticleData>  pd,
                             std::shared_ptr<ParticleGroup> pg,
                             InputFile&                     in):Base(sys,pd,pg,in){
                
                //Add bonds
                BondType::Parameters bondParameters;

                std::shared_ptr<BondType> bH_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, bH_PBC,
                                                             interactorBondParameters);

                //Constrained
                
                //ConstraintBaseType::Parameters constraintBaseParameters;
                //
                //std::shared_ptr<ConstraintBaseType> cnstr = std::make_shared<ConstraintBaseType>(this->pd,
                //                                                                                 constraintBaseParameters);
                //
                //typename ConstraintType::Parameters constraintParameters;
                //
                //constraintParameters.bondName = "BONDS";

                //constraint = std::make_shared<ConstraintType>(this->sys, this->pd, this->pg,
                //                                              this->top, cnstr,
                //                                              constraintParameters);
                
                //Add angles
                AngleType::Parameters angleParameters;

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
    
}}}}


#endif
