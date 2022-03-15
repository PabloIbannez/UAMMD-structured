#ifndef __ENM__
#define __ENM__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace ElasticNetworkModel{

    template<class Base_ >
    class ElasticNetworkModel : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::HarmonicConst_K;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            
            std::shared_ptr<InteractorBondType>   bonds;

            real K;
        
        public:
        
            ElasticNetworkModel(std::shared_ptr<System>        sys,
                                std::shared_ptr<ParticleData>  pd,
                                std::shared_ptr<ParticleGroup> pg,
                                InputFile&                     in):Base(sys,pd,pg,in),
                                                                   K(std::stof(in.getOption("ElasticNetworkModelK",InputFile::Required).str())){
                
                BondType::Parameters bondParameters;

                bondParameters.K   = K;
                
                std::shared_ptr<BondType> b = std::make_shared<BondType>(this->pd,
                                                                         bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "ENM_BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, b,
                                                             interactorBondParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
            }
    
    };
    
}}}}


#endif
