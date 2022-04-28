#ifndef __SPHERICAL_BOUND__
#define __SPHERICAL_BOUND__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Bounds{

    template<class Base_>
    class SphericalShell : public Base_{

            using Base = Base_;

        public:
            
        protected:

            using SphericalShellType           = Potentials::Bounds::SphericalShell;
            using InteractorSphericalShellType = ExternalForces<SphericalShellType>;
            
            std::shared_ptr<SphericalShellType>  sphericalPotential;
            std::shared_ptr<InteractorSphericalShellType> spherical;

        private:

            real epsilonShell;
            real sigmaShell;
            real Ashell;
            real Bshell;

            real3 shellCenter = {0,0,0};
            real  shellRadius;

        public:
        
            SphericalShell(std::shared_ptr<System>        sys,
                           std::shared_ptr<ParticleData>  pd,
                           std::shared_ptr<ParticleGroup> pg,
                           InputFile&                     in):Base(sys,pd,pg,in){

                in.getOption("epsilonShell",InputFile::Required)>>epsilonShell;
                in.getOption("sigmaShell",InputFile::Required)>>sigmaShell;
                
                in.getOption("Ashell",InputFile::Required)>>Ashell;
                in.getOption("Bshell",InputFile::Required)>>Bshell;
                
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
                                                         "Parameter Ashell added: %f",Ashell);    
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Parameter Bshell added: %f",Bshell);    
                
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Shell center: %f %f %f",
                                                         shellCenter.x,
                                                         shellCenter.y,
                                                         shellCenter.z);
                
                this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                         "Parameter shellRadius added: %f",shellRadius);    
                
                SphericalShellType::Parameters sphericalShellParameters;

                sphericalShellParameters.epsilonShell = epsilonShell;
                sphericalShellParameters.sigmaShell   = sigmaShell;
                sphericalShellParameters.Ashell       = Ashell;
                sphericalShellParameters.Bshell       = Bshell;
                
                sphericalShellParameters.shellCenter = shellCenter;
                sphericalShellParameters.shellRadius = shellRadius;
                
                {
                    auto pos = pd->getPos(access::location::cpu, access::mode::read); 
                    
                    auto groupIndex  = pg->getIndexIterator(access::location::cpu);

                    real3 dr = make_real3(pos[groupIndex[0]])-shellCenter;
                    real  innerRadius = sqrt(dot(dr,dr));

                    fori(0,this->pg->getNumberParticles()){
                              dr = make_real3(pos[groupIndex[i]])-shellCenter;
                        real  r  = sqrt(dot(dr,dr));

                        if(r>innerRadius){innerRadius=r;}
                    }

                    this->sys->template log<System::MESSAGE>("[SphericalShell] "
                                                             "Computed initial inner radius: %f",innerRadius);    
                
                    if(shellRadius < innerRadius){
                        this->sys->template log<System::CRITICAL>("[SphericalShell] "
                                                                 "Initial inner radius: %f, has to be smaller than shell radius %f",innerRadius,
                                                                                                                                    shellRadius);    
                    }
                }
                
                sphericalPotential = std::make_shared<SphericalShellType>(sphericalShellParameters);

                spherical = std::make_shared<InteractorSphericalShellType>(this->pd,this->pg,
                                                                           this->sys,
                                                                           sphericalPotential);
            }
                
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                spherical->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                spherical->updateBox(box);
            }
            
            real getSurfacePosition(){
                return shellCenter.z;
            }
    };

}}}}


#endif
