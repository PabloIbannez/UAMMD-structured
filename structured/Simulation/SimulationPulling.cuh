#ifndef __SIMULATION_PULLING__
#define __SIMULATION_PULLING__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationPulling: public Simulation<ForceField_,
                                           Minimization_,
                                           Integrator_>{

        struct partId{
            int type;
            int resId;
            int chainId;
            int modelId;
        };

        bool isPartInPullingGroup(int index,std::vector<partId>& pullingGroup){

            auto pos   = this->pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
            auto res   = this->pd->getResId(uammd::access::location::cpu,uammd::access::mode::read);
            auto chain = this->pd->getChainId(uammd::access::location::cpu,uammd::access::mode::read);
            auto model = this->pd->getModelId(uammd::access::location::cpu,uammd::access::mode::read);

            int type = pos[index].w;
            int resp = res[index];
            int chp  = chain[index];
            int mdlp = model[index];

            for(auto& pId : pullingGroup){

                if(type == pId.type    and
                   resp == pId.resId   and
                   chp  == pId.chainId and
                   mdlp == pId.modelId   ){
                    return true;
                }
            }

            return false;

        }

        void loadPullingGroup(std::shared_ptr<InputOutput::InputBlocksFile> pullingGroupsFile,
                              std::string label,std::vector<partId>& pullingGroup){
            
            //Auxiliar map which translate from name(string) to typeId(int)(internal representation)
            std::map<std::string,int> nameIdMap;

            for(int typeId : this->top->getTypes()->getTypeIdList()){
                nameIdMap[this->top->getTypes()->getTypeParameters(typeId).name]=typeId;
            }
            
            auto definitionBlock = pullingGroupsFile->getFileBlockIterator(label);
            
            std::string line;
            std::stringstream parser;

            while (definitionBlock.next(line)){

                parser.clear();
                parser.str(line);

                std::string name;  
                int resIdBuffer;
                int chnIdBuffer;
                int molIdBuffer;
                
                parser >> name         >>
                          resIdBuffer  >>
                          chnIdBuffer  >>
                          molIdBuffer  ;

                if(nameIdMap.count(name) == 0){
                    this->sys->template log<System::CRITICAL>("[SimulationPulling] Error while loading pulling group, particle type %s has not been added.", name.c_str());
                } 

                if(parser.fail()){
                    this->sys->template log<System::CRITICAL>("[SimulationPulling] Error processing line \"%s\", while reading pulling group %s.", line.c_str(),label.c_str());
                }
                
                partId partIdBuffer;

                partIdBuffer.type    = nameIdMap[name];
                partIdBuffer.resId   = resIdBuffer;
                partIdBuffer.chainId = chnIdBuffer;
                partIdBuffer.modelId = molIdBuffer;

                pullingGroup.push_back(partIdBuffer);
                
                this->sys->template log<System::DEBUG>("[SimulationPulling] Particle added to group %s, type: %s (%i), resId %i, chainId %i, modelId %i", 
                                                         label.c_str(), name.c_str(),partIdBuffer.type, partIdBuffer.resId, partIdBuffer.chainId, partIdBuffer.modelId);

            }
                
            if(pullingGroup.size() == 0){
                this->sys->template log<System::CRITICAL>("[SimulationPulling] Error while loading %s, the group is empty", label.c_str());
            } else {
                this->sys->template log<System::MESSAGE>("[SimulationPulling] Group %s loaded, size: %i", label.c_str(), pullingGroup.size());
            }
        }
    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using PullingModel = typename Interactor::ConstantForceCOMCopies;

    protected:

        std::shared_ptr<PullingModel> pullingInteractor;
        
        thrust::host_vector<int> set1;
        thrust::host_vector<int> set2;

        std::vector<std::shared_ptr<ParticleGroup>> simulationGroups;
        
        real dt;

        real pullingForce;
        
        real currentPullingForce;
        real pullingForceIncreaseRate;

        int nStepsPullingMeasure;

        std::string pullingGroupsDefinitionFilePath;

        std::string outputPullingMeasureFilePath;
        std::ofstream outputPullingMeasureFile;

    public:

        SimulationPulling(std::shared_ptr<System> sys,
                          uammd::InputFile& in,
                          bool init = true):Simulation(sys,in,false){
                                           

            in.getOption("dt",InputFile::Required)
                          >>dt;

            in.getOption("pullingForce",InputFile::Required)
                          >>pullingForce;
            in.getOption("pullingForceIncreaseRate",InputFile::Required)
                          >>pullingForceIncreaseRate;
            
            in.getOption("pullingGroupsDefinitionFilePath",InputFile::Required)
                          >>pullingGroupsDefinitionFilePath;
            
            in.getOption("nStepsPullingMeasure",InputFile::Required)
                          >>nStepsPullingMeasure;
            in.getOption("outputPullingMeasureFilePath",InputFile::Required)
                          >>outputPullingMeasureFilePath;
            
            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter dt %f"
                                                      ,dt);

            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter pullingForce %f"
                                                      ,pullingForce);
            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter pullingForceIncreaseRate %f"
                                                      ,pullingForceIncreaseRate);
            
            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter pullingGroupsDefinitionFilePath %s",
                                                      pullingGroupsDefinitionFilePath.c_str());

            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter nStepsPullingMeasure %f",
                                                      nStepsPullingMeasure);
            this->sys->template log<System::MESSAGE>("[SimulationPulling] "
                                                      "Parameter outputPullingMeasureFilePath %s",
                                                      outputPullingMeasureFilePath.c_str());
                                        
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
            
            outputPullingMeasureFile = std::ofstream(outputPullingMeasureFilePath);
            
            this->minimization->addInteractor(this->ff);
            this->integrator->addInteractor(this->ff);
            
            this->addSimulationStep(std::make_shared<InfoStep>(this->sys,this->pd,this->pg,this->nStepsInfoInterval,this->nSteps));
            this->addSimulationStep(std::make_shared<SortStep>(this->sys,this->pd,this->pg,this->nStepsSortInterval));
            
            //Pulling

            std::vector<partId> pullingGroup1;
            std::vector<partId> pullingGroup2;

            std::shared_ptr<InputOutput::InputBlocksFile> pullingGroupsFile = std::make_shared<InputOutput::InputBlocksFile>(this->sys,pullingGroupsDefinitionFilePath);

            this->loadPullingGroup(pullingGroupsFile,"PullingGroup1",pullingGroup1);
            this->loadPullingGroup(pullingGroupsFile,"PullingGroup2",pullingGroup2);

            auto id   = this->pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
            
            auto groupIndex  = this->pg->getIndexIterator(uammd::access::location::cpu);

            for(int i=0;i<this->pg->getNumberParticles();i++){

                int index = groupIndex[i];

                auto id  = this->pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
                int idp  = id[index];

                if(isPartInPullingGroup(index,pullingGroup1)){
                    set1.push_back(idp);
                }

                if(isPartInPullingGroup(index,pullingGroup2)){
                    set2.push_back(idp);
                }

            }
    
            this->sys->template log<uammd::System::MESSAGE>("[SimulationPulling] PullingSet1 size: %i, PullingSet2 size: %i",
                                                             set1.size(), set2.size());
    
            std::set<int> simIdList;
            {
                auto simId = this->pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);

                fori(0,this->pd->getNumParticles()){
                    simIdList.emplace(simId[i]);
                }
            }

            for(const int& sId : simIdList){
                selectors::simulationId selector(sId);
                simulationGroups.push_back(std::make_shared<uammd::ParticleGroup>(selector,this->pd,this->sys,"simulationId_"+std::to_string(sId)));
            }
            
            pullingInteractor = std::make_shared<PullingModel>(this->sys,
                                                               this->pd,
                                                               this->pg,
                                                               pullingGroup1.size(),pullingGroup2.size(),
                                                               set1,set2,
                                                               simIdList.size(),
                                                               PullingModel::Parameters());

            pullingInteractor->setState(real(0.0));
            
            this->integrator->addInteractor(pullingInteractor);
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

        }
        
        void next(){

            this->t++;
                
            if(abs(this->pullingInteractor->getState()) <= pullingForce*Simulation::Topology::Units::TO_INTERNAL_FORCE){
                int sign = 1;
                if(pullingForce < real(0.0)){
                    sign = -1;
                }

                this->pullingInteractor->setState(double(this->pullingInteractor->getState())+sign*double(dt*pullingForceIncreaseRate*Simulation::Topology::Units::TO_INTERNAL_FORCE));
            } else {
                this->pullingInteractor->setState(pullingForce*Simulation::Topology::Units::TO_INTERNAL_FORCE);
            }
                
            this->integrator->forwardTime();
            this->tryApplySteps();
                 
            if(this->t%this->nStepsPullingMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationPulling] Measuring indentation at step %i", this->t);
                this->measurePulling(this->outputPullingMeasureFile);
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
        
        void measurePulling(std::ofstream& out){
            out << this->t << " " << this->pullingInteractor->getState()*Simulation::Topology::Units::FROM_INTERNAL_FORCE << std::endl;
        }
};

}}

#endif
