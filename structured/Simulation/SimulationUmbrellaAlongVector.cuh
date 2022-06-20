#ifndef __SIMULATION_UMBRELLA_ALONG_VECTOR__
#define __SIMULATION_UMBRELLA_ALONG_VECTOR__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationUmbrellaAlongVector: public Simulation<ForceField_,
                                                       Minimization_,
                                                       Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using UmbrellaSetInteractor = typename Interactor::UmbrellaAlongVectorSetInteractor<>;

    protected:

        std::shared_ptr<UmbrellaSetInteractor> umbrellaSetInteractor;
        
        int Nsets;

        std::vector<std::shared_ptr<ParticleGroup>> simulationGroups;
        
        std::vector<std::shared_ptr<ParticleGroup>> simulationModelGroups1;
            
        std::vector<real3> equiPos;
        std::vector<real>  M1;

        int selectedModel;
                          
        real umbrellaK;

        real3 umbrellaInit;
        real3 umbrellaEnd;
        
        int  umbrellaWindowsNumber;
        real umbrellaWindowsSize;
        
        int umbrellaCopies;
        
        int nStepsUmbrellaMeasure;

        std::vector<std::vector<real>> rStats;

        std::string outputUmbrellaMeasureFilePath;
        std::ofstream outputUmbrellaMeasureFile;

    public:

        SimulationUmbrellaAlongVector(std::shared_ptr<System> sys,
                                      uammd::InputFile& in,
                                      bool init = true):Simulation(sys,in,false){
            
            in.getOption("umbrellaK",InputFile::Required)
                          >>umbrellaK;
            
            in.getOption("umbrellaSelectedModel",InputFile::Required)
                          >>selectedModel;
            
            in.getOption("umbrellaInit",InputFile::Required)
                          >>umbrellaInit.x >> umbrellaInit.y >> umbrellaInit.z;
            in.getOption("umbrellaEnd",InputFile::Required)
                          >>umbrellaEnd.x >> umbrellaEnd.y >> umbrellaEnd.z;
            in.getOption("umbrellaWindowsNumber",InputFile::Required)
                          >>umbrellaWindowsNumber;
            
            in.getOption("umbrellaCopies",InputFile::Required)
                          >>umbrellaCopies;
                                           
            in.getOption("nStepsUmbrellaMeasure",InputFile::Required)
                          >>nStepsUmbrellaMeasure;
            in.getOption("outputUmbrellaMeasureFilePath",InputFile::Required)
                          >>outputUmbrellaMeasureFilePath;
            
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaAlongVector] "
                                                      "Parameter umbrellaSelectedModel %i",
                                                      selectedModel);

            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaAlongVector] "
                                                      "Parameter umbrellaInit %f %f %f",
                                                      umbrellaInit.x,umbrellaInit.y,umbrellaInit.z);
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaAlongVector] "
                                                      "Parameter umbrellaEnd %f %f %f",
                                                      umbrellaEnd.x,umbrellaEnd.y,umbrellaEnd.z);
            
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaAlongVector] "
                                                      "Parameter nStepsUmbrellaMeasure %i",
                                                      nStepsUmbrellaMeasure);
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaAlongVector] "
                                                      "Parameter outputUmbrellaMeasureFilePath %s",
                                                      outputUmbrellaMeasureFilePath.c_str());
                                        
            if(init){
                this->init(in);
            }
        }

        std::vector<std::shared_ptr<ParticleGroup>> getSimulationGroups(){
            return this->simulationGroups;
        }
        
        void loadParticleBuffer(){
            this->Simulation::loadParticleBuffer();
        }
        
        void loadParticleData(){
            this->Simulation::loadParticleData();
        }

        void init(uammd::InputFile& in){
            
            //Feed particle buffer
            this->loadParticleBuffer();        
            
            //From buffer to particle data
            this->loadParticleData();
            
            //Generate group all
            this->pg = std::make_shared<ParticleGroup>(this->pd,this->sys,"All");

            //Load topology and add structure
            this->ff = std::make_shared<typename Simulation::ForceField>(this->sys,this->pd,this->pg,in);
            this->top = this->ff->getTopology();
            this->top->loadStructureData(this->pd);
            this->top->loadTypes(this->pd);
            
            this->minimization = std::make_shared<typename Simulation::Minimization>(this->pd,this->pg,this->sys,in,this->simulationStream);
            this->integrator   = std::make_shared<typename Simulation::Integrator>  (this->pd,this->pg,this->sys,in,this->simulationStream);
            
            this->integrator->template applyUnits<typename Simulation::Topology::Units>();
            
            this->minimization->addInteractor(this->ff);
            this->integrator->addInteractor(this->ff);
            
            this->addSimulationStep(std::make_shared<InfoStep>(this->sys,this->pd,this->pg,this->nStepsInfoInterval,this->nSteps));
            this->addSimulationStep(std::make_shared<SortStep>(this->sys,this->pd,this->pg,this->nStepsSortInterval));
            
            //Umbrella
            {
                
                Box box = this->integrator->getBox();

                if((box.boxSize.x < umbrellaEnd.x) or 
                   (box.boxSize.y < umbrellaEnd.y) or
                   (box.boxSize.z < umbrellaEnd.z)){
                    this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaAlongVector] Box (%f,%f,%f) has to be larger than umbrellaEnd (%f,%f,%f) in all directions.",
                                                                      box.boxSize.x,box.boxSize.y,box.boxSize.z,
                                                                      umbrellaEnd.x,umbrellaEnd.y,umbrellaEnd.z);
                }
            }
            outputUmbrellaMeasureFile = std::ofstream(outputUmbrellaMeasureFilePath);

            umbrellaWindowsSize = length(umbrellaEnd-umbrellaInit)/umbrellaWindowsNumber;
            
            this->sys->template log<uammd::System::MESSAGE>("[SimulationUmbrellaAlongVector] UmbrellaWindowsSize: %f",umbrellaWindowsSize);

            simulationGroups = groupUtils::getSimGroups(this->pd);

            int requiredSystemCopies = umbrellaWindowsNumber*umbrellaCopies;

            if(requiredSystemCopies != simulationGroups.size()){
                this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaAlongVector] The number of system copies given (%i) does not match with the number of copies required (%i).",
                                                                  simulationGroups.size(),requiredSystemCopies);
            } else {

                if(!groupUtils::checkSimGroupsEqual(this->pd,simulationGroups)){
                    this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaAlongVector] Equal simulation groups are required");
                }

                Nsets = simulationGroups.size();
                
                this->simulationModelGroups1 = groupUtils::getSimModelGroups(this->pd,selectedModel);
            }
            
            equiPos.resize(Nsets);
            M1.resize(Nsets);
            
            real3 nVector = normalize(umbrellaEnd-umbrellaInit);

            for(int i=0;i<umbrellaCopies;i++){
                for(int j=0;j<umbrellaWindowsNumber;j++){
                int setIndex = j+i*umbrellaWindowsNumber;
                equiPos[setIndex]=make_real3(j*umbrellaWindowsSize*nVector+umbrellaInit);
                
                M1[setIndex]=Measures::totalMass(simulationModelGroups1[setIndex],
                                                 this->simulationStream);
            }}

            UmbrellaSetInteractor::Parameters umbrellaParam;

            umbrellaParam.K = umbrellaK;
            
            this->umbrellaSetInteractor = std::make_shared<UmbrellaSetInteractor>(this->pg,
                                                                                  this->simulationModelGroups1,
                                                                                  equiPos,
                                                                                  umbrellaParam);
            
            this->integrator->addInteractor(this->umbrellaSetInteractor);
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
            
            //First Measure 
            this->sys->template log<System::DEBUG1>("[SimulationUmbrellaAlongVector] Measuring indentation at step %i", this->t);
            this->measureUmbrella(this->outputUmbrellaMeasureFile);

            this->t=1;

        }
        
        void next(){

            this->t++;
                
            this->integrator->forwardTime();
            this->tryApplySteps();
                 
            if(this->t%this->nStepsUmbrellaMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationUmbrellaAlongVector] Measuring indentation at step %i", this->t);
                this->measureUmbrella(this->outputUmbrellaMeasureFile);
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
        
        void measureUmbrella(std::ofstream& out){
            
            if(this->t==0){
                out << "# ";
                for(int i=0;i<umbrellaWindowsNumber;i++){
                    out << length(equiPos[i] - umbrellaInit) << " ";
                }
                out << std::endl;
                
                rStats.resize(umbrellaWindowsNumber);
                for(int i=0;i<umbrellaWindowsNumber;i++){
                    rStats[i].resize(umbrellaCopies);
                }

            }
                
            for(int i=0;i<umbrellaWindowsNumber;i++){
                rStats[i].clear();
            }
                
            Box box = this->integrator->getBox();
            for(int i=0;i<umbrellaCopies;i++){
                for(int j=0;j<umbrellaWindowsNumber;j++){
                int setIndex = j+i*umbrellaWindowsNumber;
                real m1 = M1[setIndex];
                real3 com1 = Measures::centerOfMassPos(simulationModelGroups1[setIndex],
                                                       m1,
                                                       this->simulationStream);
                real3 dr = box.apply_pbc(com1-umbrellaInit);
                rStats[j].push_back(sqrt(dot(dr,dr)));
            }}

            for(int i=0;i<umbrellaCopies;i++){
                for(int j=0;j<umbrellaWindowsNumber;j++){
                        //out << "(" << r0[j] << "," << rStats[r0[j]][i] << ") ";
                        out << rStats[j][i] << " ";
                        //std::cout << rStats[r0[j]][i] << " ";
                }
                //std::cout << std::endl;
                //std::cin.get();
                out << std::endl;
            }
        
        }
};

}}

#endif
