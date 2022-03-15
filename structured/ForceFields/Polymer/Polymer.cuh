#ifndef __POLYMER__
#define __POLYMER__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Polymer{

    template<class Base_ >
    class Polymer : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::HarmonicConst_K_r0;
            using AngleType    = Potentials::Bond3::HarmonicAngularConst_K_ang0;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;

            
            std::shared_ptr<InteractorBondType>   bonds;
            std::shared_ptr<InteractorAngleType>  angles;

        protected:

            real Kb,b0;
            real Ktheta,theta0;
        
        public:

            Polymer(std::shared_ptr<System>        sys,
                    std::shared_ptr<ParticleData>  pd,
                    std::shared_ptr<ParticleGroup> pg,
                    InputFile&                     in):Base(sys,pd,pg,in),
                                                       Kb(std::stof(in.getOption("Kb",InputFile::Required).str())),
                                                       b0(std::stof(in.getOption("b0",InputFile::Required).str())),
                                                       Ktheta(std::stof(in.getOption("Ktheta",InputFile::Required).str())),
                                                       theta0(std::stof(in.getOption("theta0",InputFile::Required).str())){
                
                this->sys->template log<System::MESSAGE>("[Polymer] Parameter Kb added: %f",Kb);
                this->sys->template log<System::MESSAGE>("[Polymer] Parameter b0 added: %f",b0);
                this->sys->template log<System::MESSAGE>("[Polymer] Parameter Ktheta added: %f",Ktheta);
                this->sys->template log<System::MESSAGE>("[Polymer] Parameter theta0 added: %f",theta0);
                
                //Add bonds
                BondType::Parameters bondParameters;

                bondParameters.K  = Kb;
                bondParameters.r0 = b0;
                
                std::shared_ptr<BondType> bH_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, bH_PBC,
                                                             interactorBondParameters);

                //Add angles
                AngleType::Parameters angleParameters;

                angleParameters.K    = Ktheta;
                angleParameters.ang0 = theta0;
                
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
