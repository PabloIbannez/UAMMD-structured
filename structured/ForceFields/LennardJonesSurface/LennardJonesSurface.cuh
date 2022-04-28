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

            using SurfaceType           = Potentials::Surface::GenericSurfacePotential;
            using InteractorSurfaceType = ExternalForces<SurfaceType>;
            
            std::shared_ptr<SurfaceType>  surfacePotential;
            std::shared_ptr<InteractorSurfaceType> surface;

        private:

            real epsilonSurf;
            real sigmaSurf;
            real surfacePosition;

        public:
        
            LennardJonesSurface(std::shared_ptr<System>        sys,
                                std::shared_ptr<ParticleData>  pd,
                                std::shared_ptr<ParticleGroup> pg,
                                InputFile&                     in):Base(sys,pd,pg,in),
                                                                    epsilonSurf(std::stof(in.getOption("epsilonSurf",InputFile::Required).str())),
                                                                    sigmaSurf(std::stof(in.getOption("sigmaSurf",InputFile::Required).str())),
                                                                    surfacePosition(std::stof(in.getOption("surfacePosition",InputFile::Required).str())){
                
                {
                    SurfaceType::Parameters surfaceParameters;

                    surfaceParameters.epsilonSurf     = epsilonSurf;
                    surfaceParameters.sigmaSurf       = sigmaSurf;
                    surfaceParameters.Asurf           = real( 1.0);
                    surfaceParameters.Bsurf           = real(-2.0);
                    surfaceParameters.offSetSurf      = real(0.0);
                    surfaceParameters.cutOffSurf      = INFINITY;
                    surfaceParameters.surfacePosition = surfacePosition;

                    surfacePotential = std::make_shared<SurfaceType>(surfaceParameters);

                    surface = std::make_shared<InteractorSurfaceType>(this->pd,this->pg,
                                                                      this->sys,
                                                                      surfacePotential);

                }
                
                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Parameter epsilonSurf added: %f",
                                                          epsilonSurf);
                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Parameter sigmaSurf added: %f",
                                                         sigmaSurf);

                this->sys->template log<System::MESSAGE>("[LennardJonesSurface] "
                                                         "Surface position added: %f",
                                                         surfacePosition);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                surface->sum(comp,st);
                Base::sum(comp,st);
            }
            
            void sum(std::string potName,Computables comp,cudaStream_t st){
                if(potName == "surface"){
                    surface->sum(comp,st);
                    return;
                }  

                this->sys->template log<System::CRITICAL>("[LennardJonesSurface] Requested potential %s to sum. "
                                                            "But %s is not present in the force field",
                                                            potName.c_str(),potName.c_str());
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
