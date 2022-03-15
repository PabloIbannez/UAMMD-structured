#ifndef __SIMULATION_UMBRELLA_COM_DISTANCE__
#define __SIMULATION_UMBRELLA_COM_DISTANCE__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationUmbrellaCenterOfMassDistance: public Simulation<ForceField_,
                                                                Minimization_,
                                                                Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using UmbrellaSetInteractor = typename Interactor::UmbrellaSetInteractor<>;

    protected:

        std::shared_ptr<UmbrellaSetInteractor> umbrellaSetInteractor;
        
        int Nsets;

        std::vector<std::shared_ptr<ParticleGroup>> simulationGroups;
        
        std::vector<std::shared_ptr<ParticleGroup>> simulationModelGroups1;
        std::vector<std::shared_ptr<ParticleGroup>> simulationModelGroups2;
            
        std::vector<real> r0;
        std::vector<real> M1;
        std::vector<real> M2;
                          
        real umbrellaK;

        real umbrellaInit;
        real umbrellaEnd;
        
        int  umbrellaWindowsNumber;
        real umbrellaWindowsSize;
        
        int umbrellaCopies;
        
        int nStepsUmbrellaMeasure;

        std::vector<std::vector<real>> rStats;

        std::string outputUmbrellaMeasureFilePath;
        std::ofstream outputUmbrellaMeasureFile;

    public:

        SimulationUmbrellaCenterOfMassDistance(std::shared_ptr<System> sys,
                                               uammd::InputFile& in,
                                               bool init = true):Simulation(sys,in,false){
            
            in.getOption("umbrellaK",InputFile::Required)
                          >>umbrellaK;
            
            in.getOption("umbrellaInit",InputFile::Required)
                          >>umbrellaInit;
            in.getOption("umbrellaEnd",InputFile::Required)
                          >>umbrellaEnd;
            in.getOption("umbrellaWindowsNumber",InputFile::Required)
                          >>umbrellaWindowsNumber;
            
            in.getOption("umbrellaCopies",InputFile::Required)
                          >>umbrellaCopies;
                                           
            in.getOption("nStepsUmbrellaMeasure",InputFile::Required)
                          >>nStepsUmbrellaMeasure;
            in.getOption("outputUmbrellaMeasureFilePath",InputFile::Required)
                          >>outputUmbrellaMeasureFilePath;
            
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaCenterOfMassDistance] "
                                                      "Parameter nStepsUmbrellaMeasure %i",
                                                      nStepsUmbrellaMeasure);
            this->sys->template log<System::MESSAGE>("[SimulationUmbrellaCenterOfMassDistance] "
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

                if((box.boxSize.x < umbrellaEnd) or 
                   (box.boxSize.y < umbrellaEnd) or
                   (box.boxSize.z < umbrellaEnd)){
                    this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaCenterOfMassDistance] Box (%f,%f,%f) has to be larger than umbrellaEnd (%f) in all directions.",
                                                                      box.boxSize.x,box.boxSize.y,box.boxSize.z,
                                                                      umbrellaEnd);
                }
            }
            outputUmbrellaMeasureFile = std::ofstream(outputUmbrellaMeasureFilePath);

            umbrellaWindowsSize = (umbrellaEnd-umbrellaInit)/umbrellaWindowsNumber;
            
            this->sys->template log<uammd::System::MESSAGE>("[SimulationUmbrellaCenterOfMassDistance] UmbrellaWindowsSize: %f",umbrellaWindowsSize);

            simulationGroups = groupUtils::getSimGroups(this->sys,this->pd);

            int requiredSystemCopies = umbrellaWindowsNumber*umbrellaCopies;

            if(requiredSystemCopies != simulationGroups.size()){
                this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaCenterOfMassDistance] The number of system copies given (%i) does not match with the number of copies required (%i).",
                                                                  simulationGroups.size(),requiredSystemCopies);
            } else {

                if(!groupUtils::checkSimGroupsEqual(this->sys,this->pd,simulationGroups)){
                    this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaCenterOfMassDistance] Equal simulation groups are required");
                }

                std::set<int> mdlList = groupUtils::getModelList(this->sys,this->pd);
                if(mdlList.size() != 2){
                    this->sys->template log<uammd::System::CRITICAL>("[SimulationUmbrellaCenterOfMassDistance] There are %i molecules for simulation but 2 are expected.",mdlList.size());
                }

                Nsets = simulationGroups.size();
                
                this->simulationModelGroups1 = groupUtils::getSimModelGroups(this->sys,this->pd,*(mdlList.begin()));
                this->simulationModelGroups2 = groupUtils::getSimModelGroups(this->sys,this->pd,*(++mdlList.begin()));
            }
            
            r0.resize(Nsets);
            M1.resize(Nsets);
            M2.resize(Nsets);
            
            for(int i=0;i<umbrellaCopies;i++){
                for(int j=0;j<umbrellaWindowsNumber;j++){
                int setIndex = j+i*umbrellaWindowsNumber;
                //r0[setIndex]=fmod(setIndex*umbrellaWindowsSize,umbrellaEnd-umbrellaInit)+umbrellaInit+umbrellaWindowsSize;
                r0[setIndex]=fmod(setIndex*umbrellaWindowsSize,umbrellaEnd-umbrellaInit)+umbrellaInit;
                
                M1[setIndex]=Measures::totalMass(this->sys,
                                                 this->pd,
                                                 simulationModelGroups1[setIndex],
                                                 this->simulationStream);
                
                M2[setIndex]=Measures::totalMass(this->sys,
                                                 this->pd,
                                                 simulationModelGroups2[setIndex],
                                                 this->simulationStream);
            }}

            UmbrellaSetInteractor::Parameters umbrellaParam;

            umbrellaParam.K = umbrellaK;
            
            this->umbrellaSetInteractor = std::make_shared<UmbrellaSetInteractor>(this->sys,
                                                                                  this->pd,
                                                                                  this->pg,
                                                                                  this->simulationModelGroups1,
                                                                                  this->simulationModelGroups2,
                                                                                  r0,
                                                                                  umbrellaParam);
            
            this->integrator->addInteractor(this->umbrellaSetInteractor);
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
            
            //First Measure 
            this->sys->template log<System::DEBUG1>("[SimulationUmbrellaCenterOfMassDistance] Measuring indentation at step %i", this->t);
            this->measureUmbrella(this->outputUmbrellaMeasureFile);

            this->t=1;

        }
        
        void next(){

            this->t++;
                
            this->integrator->forwardTime();
            this->tryApplySteps();
                 
            if(this->t%this->nStepsUmbrellaMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationUmbrellaCenterOfMassDistance] Measuring indentation at step %i", this->t);
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
                    out << r0[i] << " ";
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
                real m2 = M2[setIndex];
                real3 com1 = Measures::centerOfMassPos(this->sys,
                                                       this->pd,
                                                       simulationModelGroups1[setIndex],
                                                       m1,
                                                       this->simulationStream);
                real3 com2 = Measures::centerOfMassPos(this->sys,
                                                       this->pd,
                                                       simulationModelGroups2[setIndex],
                                                       m2,
                                                       this->simulationStream);
                //std::cout << r0[setIndex] << " set: " << setIndex << " " << com2 << " " << com1 << std::endl;
                real3 dr = box.apply_pbc(com2-com1);
                //std::cout << sqrt(dot(dr,dr)) << std::endl;
                rStats[j].push_back(sqrt(dot(dr,dr)));
            }}
            //std::cin.get();

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
