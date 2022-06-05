#ifndef __LENNARD_JONES_SURFACE__
#define __LENNARD_JONES_SURFACE__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace LennardJonesSurface{

    template<class Base_>
    class LennardJonesSurface : public Base_{

            using Base = Base_;

        public:
            
        protected:

            using SurfaceType           = Potentials::Surface::LennardJonesSurface<typename Base::Topology>;
            using InteractorSurfaceType = ExternalForces<SurfaceType>;
            
            std::shared_ptr<SurfaceType>  surfacePotential;
            std::shared_ptr<InteractorSurfaceType> surface;

        private:

            real epsilon;
            real sigma;
            real surfacePosition;

        public:
        
            LennardJonesSurface(std::shared_ptr<System>        sys,
                                std::shared_ptr<ParticleData>  pd,
                                std::shared_ptr<ParticleGroup> pg,
                                InputFile&                     in):Base(sys,pd,pg,in),
                                                                    epsilon(std::stof(in.getOption("epsilonSurf",InputFile::Required).str())),
                                                                    sigma(std::stof(in.getOption("sigmaSurf",InputFile::Required).str())),
                                                                    surfacePosition(std::stof(in.getOption("surfacePosition",InputFile::Required).str())){
                
                {
                    typename SurfaceType::Parameters surfaceParameters;

                    surfaceParameters.epsilon     = epsilon;
                    surfaceParameters.sigma       = sigma;
                    surfaceParameters.surfacePosition = surfacePosition;

                    surfacePotential = std::make_shared<SurfaceType>(surfaceParameters);

                    surface = std::make_shared<InteractorSurfaceType>(this->pd,this->pg,
                                                                      this->sys,
                                                                      surfacePotential);

                }
                
                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Parameter epsilon added: %f",
                                                          epsilon);
                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Parameter sigma added: %f",
                                                         sigma);

                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Surface position added: %f",
                                                         surfacePosition);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                surface->sum(comp,st);
            }
            
            void sum(std::string component,Computables comp,cudaStream_t st){
                if(component == "surface"){
                    surface->sum(comp,st);
                    return;
                }  
                Base::sum(component,comp,st);

                this->sys->template log<System::CRITICAL>("[LennardJonesSurface] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            component.c_str(),component.c_str());
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
            }

            real getSurfacePosition(){
                return surfacePosition;
            }
            
    };

}}}}


#endif
