#ifndef __SIMULATION_QCM__
#define __SIMULATION_QCM__

namespace uammd{
namespace structured{

template<typename ForceField_,
         typename Minimization_,
         typename Integrator_>
class SimulationQCM: public Simulation<ForceField_,
                                       Minimization_,
                                       Integrator_>{

    public:
        
        using Simulation = Simulation<ForceField_,
                                      Minimization_,
                                      Integrator_>;
        
        using QCMForcing = typename structured::Interactor::QCM_ns::QCMForcing;
        using QCM        = typename structured::Interactor::QCMInteractor;
            
    protected:

        using FixedType = Potentials::Bond1::HarmonicConst_K_r0;
        using ENMType   = Potentials::Bond2::HarmonicConst_K;
        
        using InteractorFixedType   = Interactor::BondedInteractor<FixedType,
                                                                   Interactor::BondedInteractor_ns::BondProcessor<FixedType>,
                                                                   Interactor::BondedInteractor_ns::BondReaderFromVector<FixedType>>;
        
        using InteractorENMType   = Interactor::BondedInteractor<ENMType,
                                                                 Interactor::BondedInteractor_ns::BondProcessor<ENMType>,
                                                                 Interactor::BondedInteractor_ns::BondReaderFromVector<ENMType>>;
        
        std::shared_ptr<QCMForcing> qcmf;
        std::shared_ptr<QCM>        qcm;
        
        std::shared_ptr<ParticleGroup> pgSample;

        int  qcmTypeId;
        std::vector<int> QCMids;
        std::shared_ptr<ParticleGroup> pgQCM;
        std::shared_ptr<std::vector<FixedType::Bond>> QCMfixedBonds;
        std::shared_ptr<std::vector<ENMType::Bond>>   QCMelasticNetworkModelBonds;

        std::shared_ptr<InteractorFixedType> QCMfixed;
        std::shared_ptr<InteractorENMType>   QCMenm;

        int wallTypeId;
        std::vector<int> wallIds;
        std::shared_ptr<ParticleGroup> pgWall;
        std::shared_ptr<std::vector<FixedType::Bond>> WallFixedBonds;
        std::shared_ptr<std::vector<ENMType::Bond>>   WallElasticNetworkModelBonds;
        
        std::shared_ptr<InteractorFixedType> WallFixed;
        std::shared_ptr<InteractorENMType>   WallENM;

        //Input

        Box box;

        real QCMwallMinSeparation;
        real QCMwallLatticeFactor;

        real Mqcm;
        real Rqcm;
        
        real Mwall;
        real Rwall;

        real Kqcm;
        real Kwall;
        
        /*
        real initialQCMSampleDst;
            
        real frictionConstantQCM;

        real dt;
        real descentVelocity;

        int nStepsIndentMeasure;

        std::string outputIndentationMeasureFilePath;
        
        std::ofstream outputIndentationMeasureFile;
        */

    public:

        SimulationQCM(std::shared_ptr<System> sys,
                      uammd::InputFile& in,
                      bool init = true):Simulation(sys,in,false){
                
            real3 boxSizeBuffer;

            in.getOption("boxSize",InputFile::Required)>>boxSizeBuffer.x
                                                       >>boxSizeBuffer.y
                                                       >>boxSizeBuffer.z;

            box = Box(boxSizeBuffer);
            
            in.getOption("QCMwallMinSeparation",InputFile::Required)
                          >>QCMwallMinSeparation;

            in.getOption("QCMwallLatticeFactor",InputFile::Required)
                          >>QCMwallLatticeFactor;
                    
            in.getOption("Mqcm",InputFile::Required)
                          >>Mqcm;
            in.getOption("Rqcm",InputFile::Required)
                          >>Rqcm;
            
            in.getOption("Mwall",InputFile::Required)
                          >>Mwall;
            in.getOption("Rwall",InputFile::Required)
                          >>Rwall;
 
            in.getOption("Kqcm",InputFile::Required)
                          >>Kqcm;
            in.getOption("Kwall",InputFile::Required)
                          >>Kwall;
        
            QCMfixedBonds                = std::make_shared<std::vector<FixedType::Bond>>(); 
            QCMelasticNetworkModelBonds  = std::make_shared<std::vector<ENMType::Bond>>();  
            
            WallFixedBonds               = std::make_shared<std::vector<FixedType::Bond>>();
            WallElasticNetworkModelBonds = std::make_shared<std::vector<ENMType::Bond>>();  
            
            /*
            in.getOption("frictionConstantTip",InputFile::Required)
                          >>frictionConstantTip;
                
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
             */                           
            if(init){
                this->init(in);
            }
        }
        
        void loadParticleBuffer(){
            
            this->Simulation::loadParticleBuffer();
            
            //Add QCM particles

            //fcc lattice
            real Lz = QCMwallMinSeparation*QCMwallLatticeFactor;

            int3 n;
            
            std::vector<real3> latticeVectors(3,{0,0,0});
            std::vector<real3>   basisVectors(4,{0,0,0});
            
            real3 stretchFactors;

            //Init lattice vectors
            latticeVectors[0]={1,0,0};
            latticeVectors[1]={0,1,0};
            latticeVectors[2]={0,0,1};

            //Init basis vectors
            real cellSize = 2.0*QCMwallMinSeparation/sqrt(2.0);
            basisVectors[1]={0.5,0.5,0  };
            basisVectors[2]={0.5,0  ,0.5};
            basisVectors[3]={0  ,0.5,0.5};

            //
            n.x=(int) std::ceil(box.boxSize.x/(cellSize*latticeVectors[0].x));
            n.y=(int) std::ceil(box.boxSize.y/(cellSize*latticeVectors[1].y));
            n.z=(int) std::ceil(Lz/(cellSize*latticeVectors[2].z));

            //Stretch vectors
            stretchFactors.x=box.boxSize.x/(n.x*latticeVectors[0].x);
            stretchFactors.y=box.boxSize.y/(n.y*latticeVectors[1].y);
            stretchFactors.z=Lz/(n.z*latticeVectors[2].z);

            int QCMindexStart = this->pdBuffer.size();
            for(uint i=0; i<n.x; i++){
            for(uint j=0; j<n.y; j++){
            for(uint k=0; k<n.z; k++){
                for(real3 b : basisVectors){
                    
                    int id = 0;
                    if(this->pdBuffer.size()!=0){
                        id = this->pdBuffer.back().id+1;
                    } 
                    
                    real3  pos = {-real(0.5)*box.boxSize.x + stretchFactors.x*((i + b.x)*latticeVectors[0].x+
                                                                               (j + b.y)*latticeVectors[1].x+
                                                                               (k + b.z)*latticeVectors[2].x),
                                  -real(0.5)*box.boxSize.y + stretchFactors.y*((i + b.x)*latticeVectors[0].y+
                                                                               (j + b.y)*latticeVectors[1].y+
                                                                               (k + b.z)*latticeVectors[2].y),
                                  -real(0.5)*Lz + stretchFactors.z*((i + b.x)*latticeVectors[0].z+
                                                                    (j + b.y)*latticeVectors[1].z+
                                                                    (k + b.z)*latticeVectors[2].z)};

                    this->pdBuffer.push_back({id,pos,{0,0,0}});

                    QCMids.push_back(id);
                    
                    QCMfixedBonds->push_back({id,pos});

                }
            }}}

            for(uint i=QCMindexStart; i<this->pdBuffer.size(); i++){
            for(uint j=i+1          ; j<this->pdBuffer.size(); j++){

                int id_i = this->pdBuffer[i].id;
                int id_j = this->pdBuffer[j].id;

                real3 pi = this->pdBuffer[i].pos;
                real3 pj = this->pdBuffer[j].pos;

                real3 dr = box.apply_pbc(pj-pi);
                real r0 = sqrt(dot(dr,dr));

                if(r0<QCMwallMinSeparation*real(1.1)){
                    QCMelasticNetworkModelBonds->push_back({id_i,id_j,r0});
                }
            }}
            
            int WallIndexStart = this->pdBuffer.size();
            for(uint i=0; i<n.x; i++){
            for(uint j=0; j<n.y; j++){
            for(uint k=0; k<n.z; k++){
                for(real3 b : basisVectors){
                    int    id  = this->pdBuffer.back().id+1;
                    real3  pos = {-real(0.5)*box.boxSize.x + stretchFactors.x*((i + b.x)*latticeVectors[0].x+
                                                                               (j + b.y)*latticeVectors[1].x+
                                                                               (k + b.z)*latticeVectors[2].x),
                                  -real(0.5)*box.boxSize.y + stretchFactors.y*((i + b.x)*latticeVectors[0].y+
                                                                               (j + b.y)*latticeVectors[1].y+
                                                                               (k + b.z)*latticeVectors[2].y),
                                  -real(0.5)*Lz + stretchFactors.z*((i + b.x)*latticeVectors[0].z+
                                                                    (j + b.y)*latticeVectors[1].z+
                                                                    (k + b.z)*latticeVectors[2].z)+real(0.5)*box.boxSize.z};
                    this->pdBuffer.push_back({id,pos,{0,0,0}});
                    
                    wallIds.push_back(id);
                    
                    WallFixedBonds->push_back({id,pos});

                }
            }}}
            
            for(uint i=WallIndexStart; i<this->pdBuffer.size(); i++){
            for(uint j=i+1           ; j<this->pdBuffer.size(); j++){

                int id_i = this->pdBuffer[i].id;
                int id_j = this->pdBuffer[j].id;

                real3 pi = this->pdBuffer[i].pos;
                real3 pj = this->pdBuffer[j].pos;

                real3 dr = box.apply_pbc(pj-pi);
                real r0 = sqrt(dot(dr,dr));

                if(r0<QCMwallMinSeparation*real(1.1)){
                    WallElasticNetworkModelBonds->push_back({id_i,id_j,r0});
                }
            }}
            
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
            
            //Generate sample group
            selectors::idMax sampleSel = selectors::idMax(*std::min_element(QCMids.begin(),QCMids.end()));
            this->pgSample  = std::make_shared<ParticleGroup>(sampleSel,
                                                              this->pd, this->sys, "Sample");
            
            //Generate QCM group
            selectors::idSet QCMSel = selectors::idSet(QCMids);
            this->pgQCM  = std::make_shared<ParticleGroup>(QCMSel,
                                                           this->pd, this->sys, "QCM");
            
            //Generate wall group
            selectors::idSet WallSel = selectors::idSet(wallIds);
            this->pgWall  = std::make_shared<ParticleGroup>(WallSel,
                                                           this->pd, this->sys, "Wall");
            
            //Load topology and add structure
            this-> ff = std::make_shared<typename Simulation::ForceField>(this->sys,
                                                                          this->pd,
                                                                          this->pgSample,in);
            this->top = this->ff->getTopology();
            
            {
                auto types = this->top->getTypes();

                std::vector<int> typeIdList = types->getTypeIdList();
                qcmTypeId = *(std::max_element(typeIdList.begin(),typeIdList.end()))+int(1);

                typename decltype(types)::element_type::InputTypeParameters typeBuffer;
        
                typeBuffer.name  = "QCM";
                typeBuffer.mass  =  Mqcm;
                typeBuffer.radius = Rqcm;

                types->add(qcmTypeId,typeBuffer);
                
                wallTypeId = qcmTypeId+int(1);

                typeBuffer.name  = "WALL";
                typeBuffer.mass  =  Mwall;
                typeBuffer.radius = Rwall;

                types->add(wallTypeId,typeBuffer);

                this->sys->template log<System::MESSAGE>("[SimulationQCM] QCM type Id %i",qcmTypeId);
                this->sys->template log<System::MESSAGE>("[SimulationQCM] Wall type Id %i",wallTypeId);
            }
            
            this->top->loadStructureData(this->pd);

            //QCM structure data
            {
                auto pos    = this->pd->getPos(access::location::cpu,     access::mode::write);
                auto resId  = this->pd->getResId(access::location::cpu,   access::mode::write);
                auto chnId  = this->pd->getChainId(access::location::cpu, access::mode::write);
                auto molId  = this->pd->getModelId(access::location::cpu, access::mode::write);
                
                const int *sortedIndex = this->pd->getIdOrderedIndices(access::location::cpu);
             
                for(auto id : QCMids){
                    pos[sortedIndex[id]].w = qcmTypeId;
                    resId[sortedIndex[id]] = -1;
                    chnId[sortedIndex[id]] = -1;
                    molId[sortedIndex[id]] = -1;
                }
                
                for(auto id : wallIds){
                    pos[sortedIndex[id]].w = wallTypeId;
                    resId[sortedIndex[id]] = -1;
                    chnId[sortedIndex[id]] = -1;
                    molId[sortedIndex[id]] = -1;
                }
            }

            this->top->loadTypes(this->pd);
        
            ////

            FixedType::Parameters qcmfParam;

            qcmfParam.K  = {0,Kqcm,Kqcm}; 
            qcmfParam.r0 = {0,0,0};

            std::shared_ptr<FixedType> qcmf = std::make_shared<FixedType>(this->pd,
                                                                          qcmfParam);

            InteractorFixedType::Parameters intqcmfParam;

            intqcmfParam.bondName  = "QCM_FIXED";
            this->QCMfixed = std::make_shared<InteractorFixedType>(this->sys,this->pd,
                                                                   QCMfixedBonds,qcmf,
                                                                   intqcmfParam);

            //
            
            FixedType::Parameters wallfParam;

            wallfParam.K  = {Kwall,Kwall,Kwall}; 
            wallfParam.r0 = {0,0,0};

            std::shared_ptr<FixedType> wallf = std::make_shared<FixedType>(this->pd,
                                                                           wallfParam);

            InteractorFixedType::Parameters intwallfParam;

            intwallfParam.bondName  = "WALL_FIXED";
            this->WallFixed = std::make_shared<InteractorFixedType>(this->sys,this->pd,
                                                                    WallFixedBonds,wallf,
                                                                    intwallfParam);

            ////
            
            ENMType::Parameters qcmenmParam;

            qcmenmParam.K  = Kqcm; 

            std::shared_ptr<ENMType> qcmenm = std::make_shared<ENMType>(this->pd,
                                                                        qcmenmParam);

            InteractorENMType::Parameters intqcmenmParam;

            intqcmenmParam.bondName  = "QCM_ENM";
            this->QCMenm = std::make_shared<InteractorENMType>(this->sys,this->pd,
                                                               QCMelasticNetworkModelBonds,qcmenm,
                                                               intqcmenmParam);

            //
            
            ENMType::Parameters wallenmParam;

            wallenmParam.K  = Kwall; 

            std::shared_ptr<ENMType> wallenm = std::make_shared<ENMType>(this->pd,
                                                                         wallenmParam);

            InteractorENMType::Parameters intwallenmParam;

            intwallenmParam.bondName  = "WALL_ENM";
            this->WallENM = std::make_shared<InteractorENMType>(this->sys,this->pd,
                                                                WallElasticNetworkModelBonds,wallenm,
                                                                intwallenmParam);

            ////

            this->qcmf = std::make_shared<QCMForcing>(this->sys,
                                                      in);

            this->qcmf->template applyUnits<typename Simulation::Topology::Units>();
            
            this->qcm = std::make_shared<QCM>(this->pd,
                                              this->pgQCM,
                                              this->sys,
                                              this->qcmf);
            
            this->minimization = std::make_shared<typename Simulation::Minimization>(this->pd,this->pg,this->sys,in,this->simulationStream);
            this->integrator   = std::make_shared<typename Simulation::Integrator>  (this->pd,this->pg,this->sys,in,this->simulationStream);
            
            this->integrator->template applyUnits<typename Simulation::Topology::Units>();
            
            //outputIndentationMeasureFile = std::ofstream(outputIndentationMeasureFilePath);
            
            if(this->pgSample->getNumberParticles()>0){
                this->minimization->addInteractor(this->ff);
                this->integrator->addInteractor(this->ff);
            }
            
            this->integrator->addInteractor(this->qcm);
            
            this->integrator->addInteractor(this->QCMfixed);
            this->integrator->addInteractor(this->WallFixed);
            
            this->integrator->addInteractor(this->QCMenm);
            
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

            this->t=1;

        }
        
        void next(){

            this->t++;
                
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
