#ifndef __SIMULATION_PLATES__
#define __SIMULATION_PLATES__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationPlates: public Simulation<ForceField_,
                                          Minimization_,
                                          Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using PlatesModel = typename structured::Interactor::Plates;

    protected:
        
        std::shared_ptr<PlatesModel> plates;
        
        //Input
        
        real initialTopPlatePos;
        real initialBottomPlatePos;

        real initialPlateSampleDst;
            
        real dt;
        real compressionVelocity;

        real minPlatesDst;

        int nStepsCompressionMeasure;

        std::string outputCompressionMeasureFilePath;
        
        std::ofstream outputCompressionMeasureFile;

    public:

        SimulationPlates(std::shared_ptr<System> sys,
                      uammd::InputFile& in,
                      bool init = true):Simulation(sys,in,false){
                                           
            in.getOption("initialPlateSampleDst",InputFile::Required)
                          >>initialPlateSampleDst;
            in.getOption("minimalPlatesDst",InputFile::Required)
                          >>minPlatesDst;

            in.getOption("dt",InputFile::Required)
                          >>dt;
            in.getOption("compressionVelocity",InputFile::Required)
                          >>compressionVelocity;
            
            in.getOption("nStepsCompressionMeasure",InputFile::Required)
                          >>nStepsCompressionMeasure;
            in.getOption("outputCompressionMeasureFilePath",InputFile::Required)
                          >>outputCompressionMeasureFilePath;
            
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter initialPlateSampleDst %f"
                                                      ,initialPlateSampleDst);
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter minimalPlatesDst %f"
                                                      ,minPlatesDst);
            
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter dt %f"
                                                      ,dt);
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter compressionVelocity %f"
                                                      ,compressionVelocity);

            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter nStepsCompressionMeasure %f",
                                                      nStepsCompressionMeasure);
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Parameter outputCompressionMeasureFilePath %s",
                                                      outputCompressionMeasureFilePath.c_str());
                                        
            if(init){
                this->init(in);
            }
        }
        
        void init(uammd::InputFile& in){
            Simulation::init(in);

            this->plates = std::make_shared<PlatesModel>(this->sys,
                                                         this->pd,
                                                         this->pg,in);
            
            //Set plates position
            
            real maxPartHeight = std::max_element(this->pdBuffer.begin(),
                                                  this->pdBuffer.end(),
                                                  [](const auto& a,const auto& b){return a.pos.z < b.pos.z; })->pos.z;

            this->plates->setTopPlatePosition(maxPartHeight+initialPlateSampleDst);
            
            real minPartHeight = std::max_element(this->pdBuffer.begin(),
                                                  this->pdBuffer.end(),
                                                  [](const auto& a,const auto& b){return a.pos.z > b.pos.z; })->pos.z;

            this->plates->setBottomPlatePosition(minPartHeight-initialPlateSampleDst);
            
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Initial plates positions, top: %f, bottom: %f",
                                                      maxPartHeight+initialPlateSampleDst,
                                                      minPartHeight-initialPlateSampleDst);
            
            real initialSeparation = maxPartHeight-minPartHeight+real(2.0)*initialPlateSampleDst;
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Initial plates separation: %f",
                                                      initialSeparation);
            
            int requiredStepsEstimation = std::ceil((initialSeparation-minPlatesDst)/(real(2.0)*dt*compressionVelocity));
            this->sys->template log<System::MESSAGE>("[SimulationPlates] "
                                                      "Estimated steps to get the minimal plates separation: %i",
                                                      requiredStepsEstimation);
            
            this->integrator->addInteractor(this->plates);
            
            this->outputCompressionMeasureFile = std::ofstream(outputCompressionMeasureFilePath);


        }
        
        void start(){
            
            //SimStep
            for(auto& s : this->simSteps){
                s->init(this->simulationStream);
                this->integrator->addUpdatable(s);
            }
            
            this->integrator->update();
            
            this->tryApplySteps();
            this->minimization->start();
            this->tryApplySteps();

            this->t=1;
        
            this->initialTopPlatePos=this->plates->getTopPlatePosition();
            this->initialBottomPlatePos=this->plates->getBottomPlatePosition();

        }
        
        void next(){

            this->t++;
            
            real platesSep = abs(this->plates->getTopPlatePosition()-this->plates->getBottomPlatePosition());

            if(platesSep > minPlatesDst){
                this->plates->setTopPlatePosition(initialTopPlatePos-this->t*dt*compressionVelocity);
                this->plates->setBottomPlatePosition(initialBottomPlatePos+this->t*dt*compressionVelocity);
            }

            this->integrator->forwardTime();
            this->tryApplySteps();
                 
            if(this->t%this->nStepsCompressionMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationPlates] Measuring compression at step %i", this->t);
                this->measureCompression(this->outputCompressionMeasureFile);
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
        
        real getCompressionForce(){
            real2 forces = this->plates->getPlatesForce(this->simulationStream);
            return (abs(forces.x)+abs(forces.y))/real(2.0);
        }
        
        
        void measureCompression(std::ofstream& out){
            
            real platesSep = abs(this->plates->getTopPlatePosition()-this->plates->getBottomPlatePosition());

            out << this->t   << " "
                << this->plates->getTopPlatePosition()    << " "
                << this->plates->getBottomPlatePosition() << " "
                << platesSep << " "
                << this->getCompressionForce()*Simulation::Topology::Units::FROM_INTERNAL_FORCE << std::endl;
        }
};

}}

#endif
