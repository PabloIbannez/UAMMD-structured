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
        
        using InteractorSphereType = typename Interactor::Bounds::StericSphericalShell;

    protected:
        
        std::shared_ptr<InteractorSphereType> sphere;
        
        //Input
        
        real dt;
        real compressionVelocity;
           
        real minimalSphereRadius = real(0.0);

        real initialSphereRadius;
            
    public:

        SimulationSphere(std::shared_ptr<System> sys,
                         uammd::InputFile& in,
                         bool init = true):Simulation(sys,in,false){
                                           
            in.getOption("dt",InputFile::Required)
                          >>dt;
            
            in.getOption("minimalSphereRadius",InputFile::Required)
                          >>minimalSphereRadius;
            in.getOption("compressionVelocity",InputFile::Required)
                          >>compressionVelocity;

            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter dt %f"
                                                      ,dt);

            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter minimalSphereRadius %f"
                                                      ,minimalSphereRadius);
            
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Parameter compressionVelocity %f"
                                                      ,compressionVelocity);
            
            if(init){
                this->init(in);
            }
        }
        
        void init(uammd::InputFile& in){
            Simulation::init(in);
            
            sphere = std::make_shared<InteractorSphereType>(this->pg,in);

            initialSphereRadius = sphere->getPotential()->getShellRadius(); 

            real3 sphereCenter        = sphere->getPotential()->getShellCenter(); 

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

            int requiredStepsEstimation = std::ceil((initialSphereRadius-minimalSphereRadius)/(dt*compressionVelocity));
            this->sys->template log<System::MESSAGE>("[SimulationSphere] "
                                                      "Estimated steps to get the minimal sphere radius: %i",
                                                      requiredStepsEstimation);
            
            this->minimization->addInteractor(this->sphere);
            this->integrator->addInteractor(this->sphere);
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
            
            if((this->sphere->getPotential()->getShellRadius() > minimalSphereRadius) and initialSphereRadius > minimalSphereRadius){
                this->sphere->getPotential()->setShellRadius(initialSphereRadius-dt*this->t*compressionVelocity);
            }
            if((this->sphere->getPotential()->getShellRadius() < minimalSphereRadius) and initialSphereRadius < minimalSphereRadius){
                this->sphere->getPotential()->setShellRadius(initialSphereRadius+dt*this->t*compressionVelocity);
            }

            this->integrator->forwardTime();
            this->tryApplySteps();
            
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
};

}}

#endif
