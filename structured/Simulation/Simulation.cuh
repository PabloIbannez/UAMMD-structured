#ifndef __SIMULATION__
#define __SIMULATION__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class Simulation{

    public:

        struct coordFormat{
            int   id;
            real3 pos;
            real3 vel;
            real4 dir;
        };
        
        using ForceField    = ForceField_;
        using Topology      = typename ForceField::Topology;
        
        using Minimization  = Minimization_;
        
        using Integrator    = Integrator_;

    protected:
        
        std::shared_ptr<System> sys;

        cudaStream_t simulationStream;
        
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg; //All
        
        std::shared_ptr<ForceField> ff;
        std::shared_ptr<Topology>  top;

        std::shared_ptr<Minimization> minimization;
        std::shared_ptr<Integrator>   integrator;
        
        //Buffer
        std::vector<coordFormat> pdBuffer;
        
        //SimulationSteps
        std::vector<std::shared_ptr<SimulationStep>> simSteps;

    protected:

        int t=0;

        std::string inputCoordPath;
        std::string inputTopologyPath;

        uint nSteps;
        int  nStepsInfoInterval=0;
        int  nStepsSortInterval=0;

    public:

        struct Parameters{
            std::string inputCoordPath;
            std::string inputTopologyPath;

            uint nSteps=0;
            int  nStepsInfoInterval=0;
            int  nStepsSortInterval=0;
        };
        
        static Parameters inputFileToParam(InputFile& in){

            Parameters param;
            
            in.getOption("inputCoordPath"   ,InputFile::Required)
                          >>param.inputCoordPath;
            in.getOption("inputTopologyPath",InputFile::Required)
                          >>param.inputTopologyPath;

            in.getOption("nSteps",InputFile::Optional)
                          >>param.nSteps;
            in.getOption("nStepsInfoInterval",InputFile::Optional)
                          >>param.nStepsInfoInterval;
            in.getOption("nStepsSortInterval",InputFile::Optional)
                          >>param.nStepsSortInterval;

            return param;
        }
        
        Simulation(std::shared_ptr<System> sys,
                   Parameters param,
                   uammd::InputFile& in,
                   bool init = true):sys(sys),
                                     inputCoordPath(param.inputCoordPath),
                                     inputTopologyPath(param.inputTopologyPath),
                                     nSteps(param.nSteps),
                                     nStepsInfoInterval(param.nStepsInfoInterval),
                                     nStepsSortInterval(param.nStepsSortInterval){
            
            sys->log<System::MESSAGE>("[Simulation] Parameter inputCoordPath added: %s",
                                        inputCoordPath.c_str());
            sys->log<System::MESSAGE>("[Simulation] Parameter inputTopologyPath added: %s",
                                        inputTopologyPath.c_str());
            
            sys->log<System::MESSAGE>("[Simulation] Parameter nSteps added: %i",
                                        nSteps);
            sys->log<System::MESSAGE>("[Simulation] Parameter nStepsInfoInterval added: %i",
                                        nStepsInfoInterval);
            sys->log<System::MESSAGE>("[Simulation] Parameter nStepsSortInterval added: %i",
                                        nStepsSortInterval); 
            
            cudaStreamCreate(&simulationStream);
            
            if(init){
                this->init(in);
            }
        }

        Simulation(std::shared_ptr<System> sys,
                   uammd::InputFile& in,
                   bool init = true):Simulation(sys,inputFileToParam(in),in,init){}

        ~Simulation(){
            cudaStreamDestroy(simulationStream);
        }

        cudaStream_t getSimulationStream(){return this->simulationStream;}
        
        uint getStep(){return t;}

        void setStep(uint nT){
            t=nT;
            integrator->setStep(t);
        }
        
        std::shared_ptr<ParticleData>  getParticleData(){return this->pd;}
        std::shared_ptr<ParticleGroup> getParticleGroup(){return this->pg;}
        
        std::shared_ptr<Integrator>   getIntegrator(){return this->integrator;}
        
        std::shared_ptr<Topology>   getTopology()  {return this->top;}
        std::shared_ptr<ForceField> getForceField(){return this->ff;}
        
        void addInteractor(std::shared_ptr<uammd::Interactor> interactor){
            integrator->addInteractor(interactor);
        }
        
        void addInteractorToMinimization(std::shared_ptr<uammd::Interactor> interactor){
            minimization->addInteractor(interactor);
        }

        void setConstraint(std::shared_ptr<Constraint::Constraint> constraint){
            integrator->setConstraint(constraint);
        }
        
        void addSimulationStep(std::shared_ptr<SimulationStep> ss){
            sys->log<System::DEBUG>("[Simulation] Adding simulation step %s...",ss->getName().c_str());
            simSteps.emplace_back(ss);
        }

        void tryApplySteps(){
            for(auto& s : simSteps){
                s->tryApplyStep(t,simulationStream);
            }
        }

        void init(uammd::InputFile& in){

            //Feed particle buffer
            this->loadParticleBuffer();        
            
            //From buffer to particle data
            this->loadParticleData();
            
            //Generate group all
            pg = std::make_shared<ParticleGroup>(pd,"All");

            //Load topology and add structure
            ff = std::make_shared<ForceField>(pg,in);
            top = ff->getTopology();
            top->loadStructureData(pd);
            top->loadTypes(pd);
            
            minimization = std::make_shared<Minimization>(pg,in,simulationStream);
            integrator   = std::make_shared<Integrator>  (pg,in,simulationStream);
            integrator->template applyUnits<typename Topology::Units>();
            
            minimization->addInteractor(ff);
            
            integrator->addInteractor(ff);
            if constexpr (has_isConstrained<ForceField>::value) {
                if(ff->isConstrained()){
                    integrator->setConstraint(ff->getConstraint());
                }
            }            
            this->addSimulationStep(std::make_shared<InfoStep>(pg,nStepsInfoInterval,nSteps));
            this->addSimulationStep(std::make_shared<SortStep>(pg,nStepsSortInterval));
            
        }
        
        void loadParticleBuffer(std::string inputFilePath){
            pdBuffer = InputOutput::loadCoordFromFile<coordFormat>(sys,inputFilePath);
        }
        
        void loadParticleBuffer(){
            loadParticleBuffer(inputCoordPath);
        }
        
        void updateParticleData(bool keepTypes=false){
            
            auto pId = this->pd->getId(access::location::cpu, access::mode::write);
            auto pos = this->pd->getPos(access::location::cpu, access::mode::write);
            auto vel = this->pd->getVel(access::location::cpu, access::mode::write);
            auto dir = this->pd->getDir(access::location::cpu, access::mode::write);

            fori(0,pdBuffer.size()){
                
                if(i==pdBuffer[i].id){
                    pId[i]   = pdBuffer[i].id;
                    pos[i].x = pdBuffer[i].pos.x;
                    pos[i].y = pdBuffer[i].pos.y;
                    pos[i].z = pdBuffer[i].pos.z;
                    vel[i].x = pdBuffer[i].vel.x;
                    vel[i].y = pdBuffer[i].vel.y;
                    vel[i].z = pdBuffer[i].vel.z;
                    dir[i].x = pdBuffer[i].dir.x;
                    dir[i].y = pdBuffer[i].dir.y;
                    dir[i].z = pdBuffer[i].dir.z;
                    dir[i].w = pdBuffer[i].dir.w;
                    
                    if(!keepTypes){pos[i].w = int(0);}

                } else {
                    sys->log<System::CRITICAL>("[Simulation] The internal id has to "
                                                "match with the given by coord file. Inte: %i, File %i",i,pdBuffer[i].id);
                }
            }
        }
        
        void loadParticleData(){

            //Load particles to particle data and create all group
            pd = std::make_shared<ParticleData>(pdBuffer.size(),sys);

            this->updateParticleData();
            
        }

        void next(){

            t++;

            integrator->forwardTime();
            this->tryApplySteps();
        }
        
        void start(){
            
            //SimStep
            for(auto& s : simSteps){
                s->tryInit(simulationStream);
                integrator->addUpdatable(s);
            }
            
            integrator->update();
            
            this->tryApplySteps();
            this->minimization->start();
            this->tryApplySteps();

            t=1;

        }

        void run(bool start=true){

            if(start){
                this->start();
            } else {
                for(auto& s : simSteps){
                    s->tryInit(simulationStream);
                    integrator->addUpdatable(s);
                }
                
                integrator->update();
            }
            
            Timer tim;
            tim.tic();
            
            while(t<=nSteps){
                this->next();
            }

            auto totalTime = tim.toc();

            sys->log<System::MESSAGE>("Mean FPS: %.2f", real(t)/totalTime);

        }
};

}}

#endif
