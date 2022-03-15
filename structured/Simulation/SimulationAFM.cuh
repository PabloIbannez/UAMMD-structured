#ifndef __SIMULATION_AFM__
#define __SIMULATION_AFM__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationAFM: public Simulation<ForceField_,
                                       Minimization_,
                                       Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using TipModel = typename structured::Interactor::SphericalTip;

    protected:
        
        std::shared_ptr<TipModel> tip;
        
        std::shared_ptr<ParticleGroup> pgSample;
        std::shared_ptr<ParticleGroup> pgTip;

        int tipTypeId;
        int tipId;

        //Input
        
        real Mtip;
        real Rtip;
            
        real initialTipSampleDst;
        real minimalChipHeight;
            
        real dt;
        real descentVelocity;

        int nStepsIndentMeasure;

        std::string outputIndentationMeasureFilePath;
        
        std::ofstream outputIndentationMeasureFile;

    public:

        SimulationAFM(std::shared_ptr<System> sys,
                      uammd::InputFile& in,
                      bool init = true):Simulation(sys,in,false){
                                           
            in.getOption("Mtip",InputFile::Required)
                          >>Mtip;
            in.getOption("Rtip",InputFile::Required)
                          >>Rtip;
            
            in.getOption("initialTipSampleDst",InputFile::Required)
                          >>initialTipSampleDst;
            
            in.getOption("dt",InputFile::Required)
                          >>dt;
            in.getOption("descentVelocity",InputFile::Required)
                          >>descentVelocity;
            
            in.getOption("nStepsIndentMeasure",InputFile::Required)
                          >>nStepsIndentMeasure;
            in.getOption("outputIndentationMeasureFilePath",InputFile::Required)
                          >>outputIndentationMeasureFilePath;
            
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter Mtip %f"
                                                      ,Mtip);
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter Rtip %f"
                                                      ,Rtip);
            
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter initialTipSampleDst %f"
                                                      ,initialTipSampleDst);
            
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter dt %f"
                                                      ,dt);
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter descentVelocity %f"
                                                      ,descentVelocity);

            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter nStepsIndentMeasure %f",
                                                      nStepsIndentMeasure);
            this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                      "Parameter outputIndentationMeasureFilePath %s",
                                                      outputIndentationMeasureFilePath.c_str());
            
            if(in.getOption("minimalChipHeight",InputFile::Optional)){
                in.getOption("minimalChipHeight",InputFile::Optional)>>minimalChipHeight;
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Parameter minimalChipHeight %f"
                                                          ,minimalChipHeight);
            } else {
                minimalChipHeight = std::numeric_limits<real>::lowest();
            }

                                        
            if(init){
                this->init(in);
            }
        }
        
        void loadParticleBuffer(){
            
            this->Simulation::loadParticleBuffer();

            //Add tip particles
            
            real maxPartHeight = real(0.0);
            tipId = 0;

            if(this->pdBuffer.size() > 0){
                maxPartHeight = std::max_element(this->pdBuffer.begin(),
                                                 this->pdBuffer.end(),
                                                 [](const auto& a,const auto& b){return a.pos.z < b.pos.z; })->pos.z;
                
                tipId = this->pdBuffer.back().id+1;
            }
            
            real3 tipPos = {0,0,maxPartHeight+Rtip+initialTipSampleDst};
            
            this->pdBuffer.push_back({tipId,tipPos,{0,0,0}});
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
            
            //Generate group sample
            selectors::notId sampleSel(tipId);
            this->pgSample  = std::make_shared<ParticleGroup>(sampleSel,
                                                              this->pd, this->sys, "Sample");

            //Generate group tip
            selectors::id tipSel(tipId);
            this->pgTip  = std::make_shared<ParticleGroup>(tipSel,
                                                           this->pd, this->sys, "Tip");

            //Load topology and add structure
            this-> ff = std::make_shared<typename Simulation::ForceField>(this->sys,
                                                                          this->pd,
                                                                          this->pgSample,in);
            this->top = this->ff->getTopology();

            {
                auto types = this->top->getTypes();

                std::vector<int> typeIdList = types->getTypeIdList();
                tipTypeId = *(std::max_element(typeIdList.begin(),typeIdList.end()))+int(1);

                typename decltype(types)::element_type::InputTypeParameters tipType;
        
                tipType.name  = "TIP";
                tipType.mass  =  Mtip;
                tipType.radius = Rtip;

                types->add(tipTypeId,tipType);

                this->sys->template log<System::MESSAGE>("[SimulationAFM] Tip type Id %i",tipTypeId);
            }

            this->top->loadStructureData(this->pd);

            //Tip structure data
            
            {
                auto pos    = this->pd->getPos(access::location::cpu,     access::mode::write);
                auto resId  = this->pd->getResId(access::location::cpu,   access::mode::write);
                auto chnId  = this->pd->getChainId(access::location::cpu, access::mode::write);
                auto molId  = this->pd->getModelId(access::location::cpu, access::mode::write);
                
                const int *sortedIndex = this->pd->getIdOrderedIndices(access::location::cpu);
                
                pos[sortedIndex[tipId]].w = tipTypeId;
                resId[sortedIndex[tipId]] = -1;
                chnId[sortedIndex[tipId]] = -1;
                molId[sortedIndex[tipId]] = -1;
            }

            this->top->loadTypes(this->pd);

            this->tip = std::make_shared<TipModel>(this->sys, 
                                                   this->pd, 
                                                   this->pgSample, 
                                                   this->pgTip, 
                                                   in);
            {
                real3 sampleCOM;
                
                {
                      sampleCOM = Measures::centerOfMassPos(this->sys,this->pd,
                                                            this->pgSample,
                                                            this->simulationStream);
                }
                
                auto pos = this->pd->getPos(access::location::cpu, access::mode::readwrite);
                const int *sortedIndex = this->pd->getIdOrderedIndices(access::location::cpu);

                pos[sortedIndex[tipId]].x = sampleCOM.x;
                pos[sortedIndex[tipId]].y = sampleCOM.y;
                
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Sample center of mass : %f,%f,%f",
                                                          sampleCOM.x,sampleCOM.y,sampleCOM.z);
                
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Tip position : %f,%f,%f",
                                                          pos[sortedIndex[tipId]].x,
                                                          pos[sortedIndex[tipId]].y,
                                                          pos[sortedIndex[tipId]].z);
                
                double3 chipPos = {pos[sortedIndex[tipId]].x,
                                   pos[sortedIndex[tipId]].y,
                                   pos[sortedIndex[tipId]].z};

                this->tip->setChipPosition(chipPos);
            }
            
            
            if(this->pdBuffer.size() > 0){
                real minPartHeight = std::min_element(this->pdBuffer.begin(),
                                                      this->pdBuffer.end(),
                                                      [](const auto& a,const auto& b){return a.pos.z < b.pos.z; })->pos.z;
            
                real surfacePosition = this->ff->getSurfacePosition();

                if(surfacePosition > minPartHeight){
                    this->sys->template log<System::CRITICAL>("[SimulationAFM] Surface position (%f) is larger than the minimal particle height (%f)",
                                                                surfacePosition,
                                                                minPartHeight);
                }
            }

            
            this->tip->setSurfacePosition(this->ff->getSurfacePosition());

            this->tip->template applyUnits<typename Simulation::Topology::Units>();
            
            this->minimization = std::make_shared<typename Simulation::Minimization>(this->pd,this->pg,this->sys,in,this->simulationStream);
            this->integrator   = std::make_shared<typename Simulation::Integrator>  (this->pd,this->pg,this->sys,in,this->simulationStream);
            
            {
                auto vel = this->pd->getVel(access::location::cpu, access::mode::readwrite);
                auto fricConst = this->pd->getFrictionConstant(access::location::cpu, access::mode::write);
                const int *sortedIndex = this->pd->getIdOrderedIndices(access::location::cpu);

                real frictionConstantTip = this->integrator->getFrictionConstant();
                
                auto types = this->top->getTypes();
                std::vector<int> typeIdList = types->getTypeIdList();

                real meanRadius=0;
                for(auto t : typeIdList){
                    meanRadius+=types->getTypeParameters(t).radius;
                }
                meanRadius/=types->getNumTypes();
                
                frictionConstantTip=(Rtip/meanRadius)*frictionConstantTip;
                
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Tip friction constant (computed): %f",
                                                          frictionConstantTip);

                vel[sortedIndex[tipId]] = make_real3(0,0,0);
                fricConst[sortedIndex[tipId]] = frictionConstantTip;

            }

            this->integrator->template applyUnits<typename Simulation::Topology::Units>();
            
            outputIndentationMeasureFile = std::ofstream(outputIndentationMeasureFilePath);
            
            if(minimalChipHeight != std::numeric_limits<real>::lowest()){
                int requiredStepsEstimation = std::ceil(((this->tip->getChipPosition()).z-minimalChipHeight)/(dt*descentVelocity));
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Estimated steps to get the minimal chip height: %i",
                                                          requiredStepsEstimation);
            }
            
            this->minimization->addInteractor(this->ff);
            this->integrator->addInteractor(this->ff);
            this->integrator->addInteractor(this->tip);
            
            this->addSimulationStep(std::make_shared<InfoStep>(this->sys,this->pd,this->pg,this->nStepsInfoInterval,this->nSteps));
            this->addSimulationStep(std::make_shared<SortStep>(this->sys,this->pd,this->pg,this->nStepsSortInterval));
            
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
            this->sys->template log<System::DEBUG1>("[SimulationAFM] Measuring indentation at step %i", this->t);
            this->measureIndentation(this->outputIndentationMeasureFile);

            this->t=1;

        }
        
        void next(){

            this->t++;

            if(double(this->tip->getChipHeight())>minimalChipHeight){
                this->tip->setChipHeight(double(this->tip->getChipHeight())-double(dt*descentVelocity));
            }
            this->integrator->forwardTime();
            this->tryApplySteps();
                 
            if(this->t%this->nStepsIndentMeasure==0){
                this->sys->template log<System::DEBUG1>("[SimulationAFM] Measuring indentation at step %i", this->t);
                this->measureIndentation(this->outputIndentationMeasureFile);
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
        
        real3 getTipForce(){
            //Set forces to 0
            this->integrator->resetForce();
            //Compute tip force
            this->integrator->sumForce();
            real3 tipForce = Measures::totalForce(this->sys,this->pd,pgTip,
                                                  Simulation::simulationStream); 
            
            //Set forces to 0
            this->integrator->resetForce();

            return tipForce;
        }
        
        real3 getSampleForce(){
            //Set forces to 0
            this->integrator->resetForce();
            //Compute sample force
            this->integrator->sumForce();
            real3 sampleForce = Measures::totalForce(this->sys,this->pd,pgSample, 
                                                     Simulation::simulationStream); 
            
            //Set forces to 0
            this->integrator->resetForce();

            return sampleForce;
        }
        
        void measureIndentation(std::ofstream& out){

            if(this->t==0){
                out << "# "
                    << std::setw(12) << std::left
                    << "step "   
                    << std::setw(12) << std::left
                    << "chipPos ("+Simulation::Topology::Units::L_EXTERNAL+")"   
                    << std::setw(12) << std::left
                    << "tipPos ("+Simulation::Topology::Units::L_EXTERNAL+") "   
                    << std::setw(12) << std::left
                    << "zChip ("+Simulation::Topology::Units::L_EXTERNAL+")  "   
                    << std::setw(12) << std::left
                    << "zTip ("+Simulation::Topology::Units::L_EXTERNAL+")   "   
                    << std::setw(12) << std::left
                    << "tipDfl ("+Simulation::Topology::Units::L_EXTERNAL+") "   
                    << std::setw(12) << std::left
                    << "tpFrc ("+Simulation::Topology::Units::F_EXTERNAL+")  "   
                    << std::setw(12) << std::left
                    << "smpFrc ("+Simulation::Topology::Units::F_EXTERNAL+") "   
                    << std::setw(12) << std::left
                    << "tpDflFrc ("+Simulation::Topology::Units::F_EXTERNAL+")" << std::endl;  
            }

            real3 chipPos = {real(tip->getChipPosition().x),
                             real(tip->getChipPosition().y),
                             real(tip->getChipPosition().z)};

            real3 tipPos  = tip->getTipPosition();
            
            real zChip = chipPos.z - this->ff->getSurfacePosition();
            real zTip  = tipPos.z  - this->ff->getSurfacePosition();
            
            real3 tipForce    = this->getTipForce()*Simulation::Topology::Units::FROM_INTERNAL_FORCE;
            real3 sampleForce = this->getSampleForce()*Simulation::Topology::Units::FROM_INTERNAL_FORCE;

            real tipDeflection      = tip->getTipDeflection();
            real tipDeflectionForce = tip->getTipDeflectionForce()*Simulation::Topology::Units::FROM_INTERNAL_FORCE;
                

            out << "  "  
                << std::setw(11) << std::left
                << this->t << " " 
                << std::setw(11) << std::left
                << chipPos.z << " " 
                << std::setw(11) << std::left
                << tipPos.z  << " " 
                << std::setw(11) << std::left
                << zChip     << " " 
                << std::setw(11) << std::left
                << zTip      << " " 
                << std::setw(11) << std::left
                << tipDeflection      << " " 
                << std::setw(11) << std::left
                << tipForce.z         << " "
                << std::setw(11) << std::left
                << sampleForce.z      << " " 
                << std::setw(11) << std::left
                << tipDeflectionForce << std::endl;
        }
};

}}

#endif
