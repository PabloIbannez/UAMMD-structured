#ifndef __SIMULATION_SPHERE__
#define __SIMULATION_SPHERE__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationSphere: public Simulation<ForceField_,
                                          Minimization_,
                                          Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using SphereType           = Potentials::Bounds::SphericalShell;
        using InteractorSphereType = ExternalForces<SphereType>;

    protected:
        
        std::shared_ptr<SphereType>  spherePotential;
        std::shared_ptr<InteractorSphereType> sphere;
        
        //Input
        
        real dt;
        real compressionVelocity;
            
        real epsilonSphere = real(1.0);
        real sigmaSphere   = real(1.0);
        real Asphere = real(1.0);
        real Bsphere = real(0.0);
        
        real minimalSphereRadius = real(0.0);
            
        real3 sphereCenter = {0,0,0};
        real  initialSphereRadius;
                              
        std::string outputPressureMeasureFilePath;
        std::ofstream outputPressureMeasureFile;
        
        int nStepsPressureMeasure;

    public:

        SimulationSphere(std::shared_ptr<System> sys,
                         uammd::InputFile& in,
                         bool init = true):Simulation(sys,in,false){
                                           
            in.getOption("dt",InputFile::Required)
                          >>dt;
            
            in.getOption("initialSphereRadius",InputFile::Required)
                          >>initialSphereRadius;
            in.getOption("minimalSphereRadius",InputFile::Required)
                          >>minimalSphereRadius;
            
            in.getOption("compressionVelocity",InputFile::Required)
                          >>compressionVelocity;

            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter dt %f"
                                                      ,dt);

            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter initialSphereRadius %f"
                                                      ,initialSphereRadius);
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter minimalSphereRadius %f"
                                                      ,minimalSphereRadius);
            
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter compressionVelocity %f"
                                                      ,compressionVelocity);
            
            in.getOption("epsilonSphere",InputFile::Optional)
                          >>epsilonSphere;
            in.getOption("sigmaSphere",InputFile::Optional)
                          >>sigmaSphere;
            in.getOption("Asphere",InputFile::Optional)
                          >>Asphere;
            in.getOption("Bsphere",InputFile::Optional)
                          >>Bsphere;
            
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter epsilonSphere %f"
                                                      ,epsilonSphere);
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter sigmaSphere %f"
                                                      ,sigmaSphere);
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter Asphere %f"
                                                      ,Asphere);
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter Bsphere %f"
                                                      ,Bsphere);
            
            if(in.getOption("nStepsPressureMeasure",InputFile::Optional)){
                in.getOption("nStepsPressureMeasure",InputFile::Required)
                              >>nStepsPressureMeasure;
                in.getOption("outputPressureMeasureFilePath",InputFile::Required)
                              >>outputPressureMeasureFilePath;
                
                this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                          "Parameter nStepsPressureMeasure %f",
                                                          nStepsPressureMeasure);
                this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                          "Parameter outputPressureMeasureFilePath %s",
                                                          outputPressureMeasureFilePath.c_str());
            }
            
            if(init){
                this->init(in);
            }
        }
        
        void init(uammd::InputFile& in){
            Simulation::init(in);

            {
                auto pos = this->pd->getPos(access::location::cpu, access::mode::read); 
                auto groupIndex  = this->pg->getIndexIterator(access::location::cpu);

                real3 dr = make_real3(pos[groupIndex[0]])-sphereCenter;
                real  computedInitialRadius = sqrt(dot(dr,dr));

                fori(0,this->pg->getNumberParticles()){
                          dr = make_real3(pos[groupIndex[i]])-sphereCenter;
                    real  r  = sqrt(dot(dr,dr));
                    computedInitialRadius = std::max(computedInitialRadius,r);
                }

                this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                         "Computed initial sphere radius: %f",computedInitialRadius);    

                if(initialSphereRadius < computedInitialRadius){
                    this->sys->template log<System::CRITICAL>("[SimulationSphere] "
                                                              "Computed initial sphere radius (%f) is smaller than initial radius (%f)",
                                                               computedInitialRadius,
                                                               initialSphereRadius);    
                }
            }

            SphereType::Parameters sphereParameters;

            sphereParameters.epsilonShell = epsilonSphere;
            sphereParameters.sigmaShell   = sigmaSphere;
            sphereParameters.Ashell       = Asphere;
            sphereParameters.Bshell       = Bsphere;
            
            sphereParameters.shellCenter = sphereCenter;
            sphereParameters.shellRadius = initialSphereRadius;
                                        
                
            spherePotential = std::make_shared<SphereType>(sphereParameters);

            sphere = std::make_shared<InteractorSphereType>(this->pg,
                                                            spherePotential);
            
            int requiredStepsEstimation = std::ceil((initialSphereRadius-minimalSphereRadius)/(dt*compressionVelocity));
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Estimated steps to get the minimal sphere radius: %i",
                                                      requiredStepsEstimation);
            
            this->minimization->addInteractor(this->sphere);
            this->integrator->addInteractor(this->sphere);
            
            outputPressureMeasureFile = std::ofstream(outputPressureMeasureFilePath);
        }
        
        void start(){
            
            //SimStep
            for(auto& s : this->simSteps){
                s->tryInit(this->simulationStream);
                this->integrator->addUpdatable(s);
            }
            
            this->integrator->update();
            
            this->tryApplySteps();
            this->minimization->start();
            this->tryApplySteps();

            this->t=1;
        }
        
        void next(){

            this->t++;
            
            if((this->spherePotential->getShellRadius() > minimalSphereRadius) and initialSphereRadius > minimalSphereRadius){
                this->spherePotential->setShellRadius(initialSphereRadius-dt*this->t*compressionVelocity);
            }
            if((this->spherePotential->getShellRadius() < minimalSphereRadius) and initialSphereRadius < minimalSphereRadius){
                this->spherePotential->setShellRadius(initialSphereRadius+dt*this->t*compressionVelocity);
            }

            this->integrator->forwardTime();
            this->tryApplySteps();
            
            if(this->t%this->nStepsPressureMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationSphere] Measuring pressure at step %i", this->t);
                this->measurePressure(this->outputPressureMeasureFile);
            }
        
        }
        
        void run(){

            this->start();
            
            Timer tim;
            tim.tic();
            
            while(this->t<=this->nSteps){
                this->next();
            }

            auto totalTime = tim.toc();

            this->sys->template log<System::MESSAGE>("Mean FPS: %.2f", real(this->t)/totalTime);

        }

        void measurePressure(std::ofstream& out){

        }
};

}}

#endif
