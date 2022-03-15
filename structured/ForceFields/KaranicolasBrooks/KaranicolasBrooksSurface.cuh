#ifndef __KARANICOLAS_BROOKS_SURFACE__
#define __KARANICOLAS_BROOKS_SURFACE__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace KaranicolasBrooks{

    template<class Base_>
    class KaranicolasBrooksSurface : public Base_{

            using Base = Base_;

        public:
            
        protected:

            using SurfaceType           = Potentials::Surface::KaranicolasBrooksSurfacePotential;
            using InteractorSurfaceType = ExternalForces<SurfaceType>;
            
            std::shared_ptr<SurfaceType>  surfacePotential;
            std::shared_ptr<InteractorSurfaceType> surface;

            static constexpr real rho = 1.0;
            
            static constexpr real sigmaSurface   = 3.0;
            static constexpr real epsilonSurface = 1.0;

        private:

            real chiSurface;
            real surfacePosition;

        public:
        
            KaranicolasBrooksSurface(std::shared_ptr<System>        sys,
                                     std::shared_ptr<ParticleData>  pd,
                                     std::shared_ptr<ParticleGroup> pg,
                                     InputFile&                     in):Base(sys,pd,pg,in),
                                                                    chiSurface(std::stof(in.getOption("KaranicolasBrooksSurfaceChi",InputFile::Required).str())),
                                                                    surfacePosition(std::stof(in.getOption("surfacePosition",InputFile::Required).str())){
                
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[KaranicolasBrooksSurface] Karanicolas Brooks surface force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }
                
                {
                    SurfaceType::Parameters surfaceParameters;

                    surfaceParameters.epsilonSurface  = epsilonSurface;
                    surfaceParameters.sigmaSurface    = sigmaSurface;
                    surfaceParameters.rho             = rho;
                    surfaceParameters.chiSurface      = chiSurface;
                    surfaceParameters.useInnerRadius  = true;
                    surfaceParameters.surfacePosition = surfacePosition;
                
                    surfacePotential = std::make_shared<SurfaceType>(surfaceParameters);

                    surface = std::make_shared<InteractorSurfaceType>(this->pd,this->pg,
                                                                      this->sys,
                                                                      surfacePotential);

                }
                
                this->sys->template log<System::MESSAGE>("[KaranicolasBrooksSurface] "
                                                         "Parameter KaranicolasBrooksSurfaceChi added: %f",
                                                          chiSurface);
                
                this->sys->template log<System::MESSAGE>("[KaranicolasBrooksSurface] "
                                                         "Surface position: %f",
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
