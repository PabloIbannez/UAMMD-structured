#ifndef __GENERIC_SURFACE_INTERACTOR__
#define __GENERIC_SURFACE_INTERACTOR__

namespace uammd{
namespace structured{ 
namespace Interactor{
namespace Surface{

    template<class GenericSurfaceModel> 
    class GenericSurface: public Interactor{

        protected:
            
            using SurfaceType           = Potentials::Surface::GenericSurfacePotential<GenericSurfaceModel>;
            using InteractorSurfaceType = ExternalForces<SurfaceType>;
            
            std::shared_ptr<SurfaceType>  surfacePotential;
            std::shared_ptr<InteractorSurfaceType> surface;

        public:
        
            GenericSurface(std::shared_ptr<ParticleGroup> pg,
                           InputFile&                     in):Interactor(pg,std::string("GenericSurface")){

                real epsilonSurface;
                real sigmaSurface;

                real surfacePosition;

                in.getOption("epsilonSurface",InputFile::Required)>>epsilonSurface;
                in.getOption("sigmaSurface",InputFile::Required)>>sigmaSurface;
                
                in.getOption("surfacePosition",InputFile::Required)>>surfacePosition;

                this->sys->template log<System::MESSAGE>("[GenericSurface] "
                                                         "Parameter epsilonSurface added: %f",
                                                          epsilonSurface);
                this->sys->template log<System::MESSAGE>("[GenericSurface] "
                                                         "Parameter sigmaSurface added: %f",
                                                         sigmaSurface);

                this->sys->template log<System::MESSAGE>("[GenericSurface] "
                                                         "Parameter surfacePosition added: %f",
                                                         surfacePosition);    
                
                typename SurfaceType::Parameters genericSurfaceParameters;

                genericSurfaceParameters.epsilon = epsilonSurface;
                genericSurfaceParameters.sigma   = sigmaSurface;
                
                genericSurfaceParameters.surfacePosition = surfacePosition;
                
                surfacePotential = std::make_shared<SurfaceType>(genericSurfaceParameters);

                surface = std::make_shared<InteractorSurfaceType>(this->pg,
                                                                  surfacePotential);
            }

            std::shared_ptr<SurfaceType> getPotential(){
                return surfacePotential;
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                surface->sum(comp,st);
            }
            
            void updateBox(Box box){
                surface->updateBox(box);
            }
    };

    using StericSurface       = GenericSurface<Potentials::Surface::GenericSurfacePotential_ns::GenericSurfacePotentialModel::Steric>; 
    using LennardJonesSurface = GenericSurface<Potentials::Surface::GenericSurfacePotential_ns::GenericSurfacePotentialModel::LennardJones>; 

}}}}


#endif
