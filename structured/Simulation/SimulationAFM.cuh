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
            
        real dt;
        real descentVelocity;
            
        bool setTargetIndentationForce;
        real targetIndentationForce;
        bool targetIndentationForceReached;
        
        bool setMinimalChipHeight;
        real minimalChipHeight;
        bool minimalChipHeightReached;
        
        bool checkSurfacePosition;

        int nStepsIndentMeasure;
        
        std::string outputIndentationMeasureFilePath;
        std::ofstream outputIndentationMeasureFile;

    public:

        struct Parameters: public Simulation::Parameters{
            real Mtip;
            real Rtip;
                
            real initialTipSampleDst;
                
            real dt;
            real descentVelocity;
            
            bool setMinimalChipHeight;
            real minimalChipHeight;
            
            bool setTargetIndentationForce;
            real targetIndentationForce;
            
            bool checkSurfacePosition=true;

            int nStepsIndentMeasure;

            std::string outputIndentationMeasureFilePath;
        };
        
        static Parameters inputFileToParam(InputFile& in){

            Parameters param;
            static_cast<typename Simulation::Parameters&>(param) = Simulation::inputFileToParam(in); 
            
            in.getOption("Mtip",InputFile::Required)
                          >>param.Mtip;
            in.getOption("Rtip",InputFile::Required)
                          >>param.Rtip;
            
            in.getOption("initialTipSampleDst",InputFile::Required)
                          >>param.initialTipSampleDst;
            
            in.getOption("dt",InputFile::Required)
                          >>param.dt;
            in.getOption("descentVelocity",InputFile::Required)
                          >>param.descentVelocity;
            
            in.getOption("nStepsIndentMeasure",InputFile::Required)
                          >>param.nStepsIndentMeasure;
            in.getOption("outputIndentationMeasureFilePath",InputFile::Required)
                          >>param.outputIndentationMeasureFilePath;
            
            if(in.getOption("minimalChipHeight",InputFile::Optional)){
                in.getOption("minimalChipHeight",InputFile::Optional)>>param.minimalChipHeight;
                param.setMinimalChipHeight = true;
            } else {
                param.minimalChipHeight    = std::numeric_limits<real>::lowest();
                param.setMinimalChipHeight = false;
            }
            
            if(in.getOption("targetIndentationForce",InputFile::Optional)){
                in.getOption("targetIndentationForce",InputFile::Optional)>>param.targetIndentationForce;
                param.setTargetIndentationForce = true;
            } else {
                param.targetIndentationForce = std::numeric_limits<real>::max();
                param.setTargetIndentationForce = false;
            }
            
            if(in.getOption("notCheckSurfacePosition",InputFile::Optional)){
                param.checkSurfacePosition = false;
            } else {
                param.checkSurfacePosition = true;
            }

            return param;
        }
        
        SimulationAFM(std::shared_ptr<System> sys,
                      Parameters param,
                      uammd::InputFile& in,
                      bool init = true):Simulation(sys,param,in,false),
                                        Mtip(param.Mtip),
                                        Rtip(param.Rtip),
                                        initialTipSampleDst(param.initialTipSampleDst),
                                        dt(param.dt),
                                        descentVelocity(param.descentVelocity),
                                        nStepsIndentMeasure(param.nStepsIndentMeasure),
                                        outputIndentationMeasureFilePath(param.outputIndentationMeasureFilePath),
                                        setTargetIndentationForce(param.setTargetIndentationForce),
                                        targetIndentationForce(param.targetIndentationForce),
                                        targetIndentationForceReached(false),
                                        setMinimalChipHeight(param.setMinimalChipHeight),
                                        minimalChipHeight(param.minimalChipHeight),
                                        minimalChipHeightReached(false),
                                        checkSurfacePosition(param.checkSurfacePosition){
                
            
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
            
            if(setMinimalChipHeight){
                this->sys->template log<System::MESSAGE>("[SimulationAFM] "
                                                          "Parameter minimalChipHeight %f"
                                                          ,minimalChipHeight);
            }
            
            if(setTargetIndentationForce){
                this->sys->template log<System::WARNING>("[SimulationAFM] "
                                                         "Target force set to:%f",targetIndentationForce);
            }
            
            if(!checkSurfacePosition){
                this->sys->template log<System::WARNING>("[SimulationAFM] "
                                                         "Not performing surface position checking...");
            }
            
            if(init){
                this->init(in);
            }

        }

        SimulationAFM(std::shared_ptr<System> sys,
                      uammd::InputFile& in,
                      bool init = true):SimulationAFM(sys,inputFileToParam(in),in,init){}
        
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
            this->pg = std::make_shared<ParticleGroup>(this->pd,"All");
            
            //Generate group sample
            selectors::notId sampleSel(tipId);
            this->pgSample  = std::make_shared<ParticleGroup>(sampleSel,
                                                              this->pd,"Sample");

            //Generate group tip
            selectors::id tipSel(tipId);
            this->pgTip  = std::make_shared<ParticleGroup>(tipSel,
                                                           this->pd,"Tip");

            //Load topology and add structure
            this-> ff = std::make_shared<typename Simulation::ForceField>(this->pgSample,in);

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

            this->tip = std::make_shared<TipModel>(this->pgSample, 
                                                   this->pgTip, 
                                                   in);
            {
                real3 sampleCOM;
                
                {
                      sampleCOM = Measures::centerOfMassPos(this->pgSample,
                                                            Measures::totalMass(this->pgSample,this->simulationStream),
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
            
            
            if(checkSurfacePosition){
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
            }
            
            this->tip->setSurfacePosition(this->ff->getSurfacePosition());

            this->tip->template applyUnits<typename Simulation::Topology::Units>();
            
            this->minimization = std::make_shared<typename Simulation::Minimization>(this->pg,in,this->simulationStream);
            this->integrator   = std::make_shared<typename Simulation::Integrator>  (this->pg,in,this->simulationStream);
            
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
            
            this->minimization->addInteractor(this->ff);
            this->integrator->addInteractor(this->ff);
            this->integrator->addInteractor(this->tip);
            
            this->addSimulationStep(std::make_shared<InfoStep>(this->pg,this->nStepsInfoInterval,this->nSteps));
            this->addSimulationStep(std::make_shared<SortStep>(this->pg,this->nStepsSortInterval));
            
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
            this->sys->template log<System::DEBUG1>("[SimulationAFM] Measuring indentation at step %i", this->t);
            this->measureIndentation(this->outputIndentationMeasureFile);

            this->t=1;

        }
        
        void next(){

            this->t++;

            this->tip->setChipHeight(double(this->tip->getChipHeight())-double(dt*descentVelocity));

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
            
            if(!setTargetIndentationForce and !setMinimalChipHeight){
                while(this->t<=this->nSteps){
                    this->next();
                }
            } else {
                while(!targetIndentationForceReached and !minimalChipHeightReached){
                    this->next();
                }
            }

            auto totalTime = tim.toc();

            this->sys->template log<System::MESSAGE>("Mean FPS: %.2f", real(this->t)/totalTime);

        }
        
        std::tuple<real3,real3> getTipAndSampleForce(){

            uninitialized_cached_vector<real4> forceBuffer(this->pg->getNumberParticles());

            //Copy force to buffer
            {
                auto force = this->pd->getForce(access::location::gpu, access::mode::read);     
                thrust::copy(thrust::cuda::par.on(Simulation::simulationStream),
                             force.begin(), 
                             force.end(), 
                             forceBuffer.begin());
            }

            //Set forces to 0
            this->integrator->resetForce();
            //Compute tip force
            this->integrator->updateForce();
            real3 tipForce = Measures::totalForce(pgTip,
                                                  Simulation::simulationStream); 
            real3 sampleForce = Measures::totalForce(pgSample, 
                                                     Simulation::simulationStream); 
            
            //Copy back buffer to force
            {
                auto force = this->pd->getForce(access::location::gpu, access::mode::read);     
                thrust::copy(thrust::cuda::par.on(Simulation::simulationStream),
                             forceBuffer.begin(), 
                             forceBuffer.end(), 
                             force.begin());
            }

            return {tipForce,sampleForce};
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
            
            auto tipAndSampleForce = this->getTipAndSampleForce();

            real3 tipForce    = std::get<0>(tipAndSampleForce)*Simulation::Topology::Units::FROM_INTERNAL_FORCE;
            real3 sampleForce = std::get<1>(tipAndSampleForce)*Simulation::Topology::Units::FROM_INTERNAL_FORCE;

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

            if(chipPos.z < minimalChipHeight){
                minimalChipHeightReached = true;
            }

            if(tipDeflectionForce > targetIndentationForce){
                targetIndentationForceReached = true;
            }
        }
};

}}

#endif
