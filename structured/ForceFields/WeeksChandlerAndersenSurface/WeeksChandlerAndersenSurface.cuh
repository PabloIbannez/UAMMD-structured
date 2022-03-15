#ifndef __WCA_SURFACE__
#define __WCA_SURFACE__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace WeeksChandlerAndersenSurface{

    template<class Base_>
    class WeeksChandlerAndersenSurface : public Base_{

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
        
            WeeksChandlerAndersenSurface(std::shared_ptr<System>        sys,
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
                    surfaceParameters.offSetSurf      = epsilonSurf;
                    surfaceParameters.cutOffSurf      = sigmaSurf;
                    surfaceParameters.surfacePosition = surfacePosition;

                    surfacePotential = std::make_shared<SurfaceType>(surfaceParameters);

                    surface = std::make_shared<InteractorSurfaceType>(this->pd,this->pg,
                                                                      this->sys,
                                                                      surfacePotential);

                }
                
                this->sys->template log<System::MESSAGE>("[WeeksChandlerAndersenSurface] "
                                                         "Parameter epsilonSurf added: %f",
                                                          epsilonSurf);
                this->sys->template log<System::MESSAGE>("[WeeksChandlerAndersenSurface] "
                                                         "Parameter sigmaSurf added: %f",
                                                         sigmaSurf);

                this->sys->template log<System::MESSAGE>("[WeeksChandlerAndersenSurface] "
                                                         "Surface position added: %f",
                                                         surfacePosition);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                surface->sum(comp,st);
                Base::sum(comp,st);
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
