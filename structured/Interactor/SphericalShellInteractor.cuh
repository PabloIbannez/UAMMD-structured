#ifndef __SPHERICAL_SHELL_INTERACTOR__
#define __SPHERICAL_SHELL_INTERACTOR__

namespace uammd{
namespace structured{ 
namespace Interactor{
namespace Bounds{

    template<class SphericalShellModel> 
    class SphericalShell : public Interactor{

        protected:
            
            using SphericalShellType           = Potentials::Bounds::SphericalShellPotential<SphericalShellModel>;
            using InteractorSphericalShellType = ExternalForces<SphericalShellType>;
            
            std::shared_ptr<SphericalShellType>  sphericalPotential;
            std::shared_ptr<InteractorSphericalShellType> spherical;

        public:
        
            SphericalShell(std::shared_ptr<ParticleGroup> pg,
                           InputFile&                     in):Interactor(pg,std::string("SphericalShape")){

                real epsilonShell;
                real sigmaShell;

                real3 shellCenter = {0,0,0};
                real  shellRadius;

                in.getOption("epsilonShell",InputFile::Required)>>epsilonShell;
                in.getOption("sigmaShell",InputFile::Required)>>sigmaShell;
                
                in.getOption("shellCenter",InputFile::Required)>>shellCenter.x>>
                                                                 shellCenter.y>>
                                                                 shellCenter.z;
                
                in.getOption("shellRadius",InputFile::Required)>>shellRadius;
                
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Parameter epsilonShell added: %f",
                                                          epsilonShell);
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Parameter sigmaShell added: %f",
                                                         sigmaShell);

                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Shell center: %f %f %f",
                                                         shellCenter.x,
                                                         shellCenter.y,
                                                         shellCenter.z);
                
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Parameter shellRadius added: %f",shellRadius);    
                
                typename SphericalShellType::Parameters sphericalShellParameters;

                sphericalShellParameters.epsilonShell = epsilonShell;
                sphericalShellParameters.sigmaShell   = sigmaShell;
                
                sphericalShellParameters.shellCenter = shellCenter;
                sphericalShellParameters.shellRadius = shellRadius;
                
                sphericalPotential = std::make_shared<SphericalShellType>(sphericalShellParameters);

                spherical = std::make_shared<InteractorSphericalShellType>(this->pg,
                                                                           sphericalPotential);
            }

            std::shared_ptr<SphericalShellType> getPotential(){
                return sphericalPotential;
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                spherical->sum(comp,st);
            }
            
            void updateBox(Box box){
                spherical->updateBox(box);
            }
    };

    using StericSphericalShell       = SphericalShell<Potentials::Bounds::SphericalShellPotential_ns::SphericalShellPotentialModel::Steric>; 
    using LennardJonesSphericalShell = SphericalShell<Potentials::Bounds::SphericalShellPotential_ns::SphericalShellPotentialModel::LennardJones>; 

}}}}


#endif
