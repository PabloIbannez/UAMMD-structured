#ifndef __CLASHED__
#define __CLASHED__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Clashed{

    template<class Intra_  >
    class Clashed : public Intra_{

            using Base = Intra_;
        
        protected:

            using CLHDType    = Potentials::UnBound::Clashed;

            using InteractorCLHDType   = Interactor::PairInteractor<CLHDType,typename Base::NeighbourList>;

            std::shared_ptr<CLHDType>         potCLHD;
            std::shared_ptr<InteractorCLHDType>  clhd;
        
        protected:

            real lambda;
            real gamma;
            
            real cutOffDstCLHD;

        public:
        
            Clashed(std::shared_ptr<System>        sys,
                    std::shared_ptr<ParticleData>  pd,
                    std::shared_ptr<ParticleGroup> pg,
                    InputFile&                     in):Base(sys,pd,pg,in),
                                                       lambda(std::stof(in.getOption("lambda",InputFile::Required).str())),
                                                       gamma(std::stof(in.getOption("gamma",InputFile::Required).str())){

                real maxRadius = real(0.0);
                
                auto typesparam = this->top->getTypes();
                for(int typeId : typesparam->getTypeIdList()){
                    real radius = typesparam->getTypeParameters(typeId).radius;
                    if(maxRadius < radius){
                        maxRadius = radius;
                    }
                }

                cutOffDstCLHD = real(2.0)*gamma*maxRadius*real(1.01);
                    
                this->sys->template log<System::MESSAGE>("[Clashed] "
                                                         "Parameter cutOffDstCLHD computed: %f",
                                                          cutOffDstCLHD);
                
                if(cutOffDstCLHD >= this->nl->getCutOffVerlet()){
                    this->sys->template log<System::CRITICAL>("[Clashed] cutOffDstCLHD (%f) "
                                                              "has to be smaller than VerletListDst (%f)",
                                                               cutOffDstCLHD,this->nl->getCutOffVerlet());
                }

                this->nl->setCutOff(cutOffDstCLHD);
                
                this->sys->template log<System::MESSAGE>("[Clashed] "
                                                         "Parameter lambda added: %f",
                                                          lambda);

                this->sys->template log<System::MESSAGE>("[Clashed] "
                                                         "Parameter gamma added: %f",
                                                          gamma);
                
                //Add clhd
                typename CLHDType::Parameters clhdPotentialParam;
                
                clhdPotentialParam.lambda = lambda;
                clhdPotentialParam.gamma  = gamma;
                
                clhdPotentialParam.cutOffDst = cutOffDstCLHD;
                
                potCLHD = std::make_shared<CLHDType>(this->pd,clhdPotentialParam);

                typename InteractorCLHDType::Parameters interactorCLHDParameters;

                interactorCLHDParameters.name = "CLHD";
                interactorCLHDParameters.pot  = potCLHD;
                interactorCLHDParameters.nl   = this->nl;
                interactorCLHDParameters.conditionInteractionName = "inter";

                clhd = std::make_shared<InteractorCLHDType>(this->sys,this->pd,this->pg,
                                                            interactorCLHDParameters);

            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                clhd->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                clhd->updateBox(box);
            }
            
    };

}}}}


#endif
