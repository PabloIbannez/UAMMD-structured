#ifndef __SIMULATION_STEP__
#define __SIMULATION_STEP__

namespace uammd{
namespace structured{

class SimulationStep: public ParameterUpdatable{
    
    protected:

        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;

        std::string name;

        bool initialized = false;

        int interval;
        int skip;

        int lastStepApplied;
        
        virtual void init( cudaStream_t st) = 0;
        virtual void applyStep(int step, cudaStream_t st) = 0;
            
    public:

        SimulationStep(std::shared_ptr<ParticleGroup> pg,
                       std::string name,
                       int interval,
                       int skip):pg(pg),
                                 pd(pg->getParticleData()),
                                 sys(pg->getParticleData()->getSystem()),
                                 name(name),
                                 interval(interval),
                                 skip(skip),
                                 lastStepApplied(-1){}
        
        SimulationStep(std::shared_ptr<ParticleGroup> pg,
                       std::string name,
                       int interval):SimulationStep(pg,name,interval,0){}

        ~SimulationStep(){}

        std::string getName(){return name;}
        int getLastStepApplied(){return lastStepApplied;}

        virtual void tryInit( cudaStream_t st){
            if(!initialized){
                this->init(st);
                initialized=true;
            } 
        };

        void tryInit(){
            cudaDeviceSynchronize();
            tryInit(0);
            cudaDeviceSynchronize();
        }

        virtual void tryApplyStep(int step, cudaStream_t st,bool force=false){
            if(initialized){
                if(this->interval==0 and !force){return;}
                if(((step%this->interval)==0 and (step>=skip)) or force){
                    lastStepApplied = step;
                    this->applyStep(step,st);
                }
            }
        }
        
        void tryApplyStep(int step,bool force=false){
            cudaDeviceSynchronize();
            tryApplyStep(step,0,force);
            cudaDeviceSynchronize();
        }

};

class SortStep: public SimulationStep {

    public:
        
    SortStep(std::shared_ptr<ParticleGroup> pg,
             int interval):SimulationStep(pg,"SortStep",interval){}

        void init(cudaStream_t st) override {}

        void applyStep(int step, cudaStream_t st) override {
            this->pd->sortParticles(st);
        }

};

class InfoStep: public SimulationStep {
    
    public:
    
    Timer tim;

    int nSteps;

    InfoStep(std::shared_ptr<ParticleGroup> pg,
             int interval,
             int nSteps):SimulationStep(pg,"InfoStep",interval),nSteps(nSteps){}

    void init(cudaStream_t st) override{
        tim = Timer();
        tim.tic();
    }

    void applyStep(int step, cudaStream_t st) override{

        if(step != 0){
            real time = tim.toc();

            if(nSteps != 0){
                std::chrono::seconds rtime_s(uint((nSteps-step)/(real(step)/time)));

                int hours   = std::chrono::duration_cast<std::chrono::hours>(rtime_s).count();
                int minutes = std::chrono::duration_cast<std::chrono::minutes>(rtime_s).count()%60;
                int seconds = std::chrono::duration_cast<std::chrono::seconds>(rtime_s).count()%60;
                
                this->sys->template log<System::MESSAGE>("Step %i, ETA: h:%i m:%i s:%i, mean FPS: %0.2f", step, hours,minutes,seconds,real(step)/time);
            } else {
                this->sys->template log<System::MESSAGE>("Step %i, mean FPS: %0.2f", step, real(step)/time);
            }

        }
    }

};

template <class Units>
class WriteStep: public SimulationStep{

    public:

        std::vector<std::string> outPutFormatAvail={"coord","sp","xyz","pdb","itpv","itpd","dcd","lammpstrj","back","vel"};
        
        struct Parameters{

            bool pbc    = true;
            bool append = false;

            std::string outPutFormat;
            std::string outPutFilePath;
        };
    
    private:
        
        std::string outPutFileName;
        std::string outPutFilePath;
        std::string outPutFormat;

        std::ofstream outPutFile;

        bool pbc;
        bool append;

        int lastWrittenStep;
        int frame;

        real dt;
        real time;
        Box box;
        
    public:
        
        static Parameters inputFileToParam(InputFile& in){

            Parameters param;

            param.pbc    = true;
            param.append = false;

            in.getOption("outPutFormat"  ,InputFile::Required)
                          >>param.outPutFormat;
            in.getOption("outPutFilePath",InputFile::Required)
                          >>param.outPutFilePath;
            
            return param;
        }
        
        WriteStep(std::shared_ptr<ParticleGroup> pg,
                  int interval,
                  bool pbc,
                  bool append,
                  std::string outPutFormat,
                  std::string outPutFilePath):WriteStep(pg,interval,{pbc,append,outPutFormat,outPutFilePath}){}
        
        WriteStep(std::shared_ptr<ParticleGroup> pg,
                  int interval,
                  InputFile& in):WriteStep(pg,interval,inputFileToParam(in)){}
        
        WriteStep(std::shared_ptr<ParticleGroup> pg,
                  InputFile& in):WriteStep(pg,
                                           std::stoi(in.getOption("nStepsWriteInterval",InputFile::Required).str()),
                                           inputFileToParam(in)){}

        WriteStep(std::shared_ptr<ParticleGroup> pg,
                  int interval,
                  Parameters param):SimulationStep(pg,"Output",interval),
                                    pbc(param.pbc),
                                    append(param.append),
                                    outPutFileName(param.outPutFilePath),
                                    outPutFormat(param.outPutFormat),
                                    outPutFilePath(param.outPutFilePath){}
    
        void init(cudaStream_t st) override{

            lastWrittenStep = -1;
            frame=0;

            sys->log<System::DEBUG>("[WriteStep] Parameter outPutFormat added: %s",
                                     outPutFormat.c_str());
            sys->log<System::DEBUG>("[WriteStep] Parameter outPutFilePath added: %s",
                                     outPutFilePath.c_str());
            
            if(this->interval != 0){
                if(std::find(outPutFormatAvail.begin(),
                             outPutFormatAvail.end(),
                             outPutFormat)==outPutFormatAvail.end()){

                    std::string message = "[WriteStep] The selected output format ("+outPutFormat+") is not available,options are:";

                    int i;
                    for(i=0;i<outPutFormatAvail.size()-1;i++){
                        message.append(" "+outPutFormatAvail[i]+",");
                    }
                    message.append(" "+outPutFormatAvail[i]+".");

                    sys->log<System::CRITICAL>("%s", message.c_str());
                } else {
                    
                    if(outPutFormat == 
                       std::string("dcd")){
                        
                        if(!append){
                            std::ofstream outPSF(outPutFilePath + ".psf");
                            
                            psf::WritePSF(pg,outPSF);

                            outPutFilePath = outPutFilePath + "." + outPutFormat;
                            outPutFile.open(outPutFilePath,std::ios::binary);
                    
                            InputOutput::Output::WriteDCDheader(pg,
                                                                frame,
                                                                this->interval,
                                                                outPutFile);
                        } else {
                            outPutFilePath = outPutFilePath + "." + outPutFormat;
                            outPutFile.open(outPutFilePath,std::ios::binary | std::ios_base::app);
                        }
                    } else {
                        outPutFilePath = outPutFilePath + "." + outPutFormat;
                        if(!append){
                            outPutFile.open(outPutFilePath);
                        } else {
                            outPutFile.open(outPutFilePath,std::ios_base::app);
                        }
                    }
                }
            } else {
                sys->log<System::WARNING>("[WriteStep] Write interval is 0, no files will be generated!");
            }
            
        }

        void applyStep(int step, cudaStream_t st) override{

            if(step==lastWrittenStep){
                return;
            } else {
                lastWrittenStep = step;
            }

            frame++;

            if(outPutFormat == 
               std::string("coord")){
                InputOutput::Output::WriteCoord(pg,box,outPutFile);
            } else if(outPutFormat == 
                      std::string("sp")){
                InputOutput::Output::WriteSP(pg,box,outPutFile);
            } else if (outPutFormat == 
                       std::string("xyz")){
                InputOutput::Output::WriteXYZ(pg,box,outPutFile);
            } else if (outPutFormat == 
                       std::string("pdb")){
                InputOutput::Output::WritePDB(pg,box,frame,outPutFile);
            } else if (outPutFormat == 
                       std::string("itpv")){
                InputOutput::Output::WriteITPV(pg,box,outPutFile);
            } else if (outPutFormat == 
                       std::string("itpd")){
                InputOutput::Output::WriteITPD(pg,box,outPutFile);
            } else if (outPutFormat == 
                       std::string("dcd")){
                InputOutput::Output::WriteDCD(pg,
                                              box,
                                              frame,step,
                                              outPutFile);
            } else if (outPutFormat == 
                       std::string("lammpstrj")){
                InputOutput::Output::WriteLAMMPS(pg,box,time,outPutFile);
            } else if (outPutFormat == 
                       std::string("vel")){
                InputOutput::Output::WriteVelocity(pg,outPutFile);
            } else if (outPutFormat == 
                       std::string("back")){
                outPutFile.close();
                outPutFile.open(outPutFilePath);
                InputOutput::Output::WriteCoord(pg,box,outPutFile);
            } else {
                sys->log<System::CRITICAL>("[WriteStep] It should not have happen,"
                                           "have the output format changed ?");
            }
        }
            
        virtual void updateBox(Box newBox) override {
            box = newBox;

            if(!pbc){
                box.setPeriodicity(false,false,false);
            }
        }
        
        virtual void updateTimeStep(real newDt) override {
            dt = newDt*Units::FROM_INTERNAL_TIME;
        }
        
        virtual void updateSimulationTime(real simulationTime) override {
            time = simulationTime*Units::FROM_INTERNAL_TIME;
        }

        void setPBC(bool new_pbc){
            pbc=new_pbc;
        }
};

template<class ForceField>
class ComputableMeasure: public SimulationStep{

        std::ofstream outPutFile;

        std::shared_ptr<ForceField> ff;

        uammd::Interactor::Computables compToMeasure;

    public:
        
        ComputableMeasure(std::shared_ptr<ParticleGroup> pg,
                          int interval,
                          std::string outPutFileName,
                          std::shared_ptr<ForceField> ff,
                          uammd::Interactor::Computables compToMeasure):SimulationStep(pg,"ComputableMeasure",interval),
                                                                        ff(ff),compToMeasure(compToMeasure){
                outPutFile = std::ofstream(outPutFileName);
        }

        void init(cudaStream_t st) override{}

        void applyStep(int step, cudaStream_t st) override{

            if(compToMeasure.energy == true){
                uninitialized_cached_vector<real> energyBuffer(this->pg->getNumberParticles());

                //Copy energy to buffer
                {
                    auto energy = this->pd->getEnergy(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 energy.begin(), 
                                 energy.end(), 
                                 energyBuffer.begin());
                }

                this->integrator->resetEnergy();
                this->integrator->updateEnergy();
                real totalEnergy = Measures::totalPotentialEnergy(pg, 
                                                                  st); 
                
                //Copy back buffer to energy
                {
                    auto energy = this->pd->getEnergy(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 energyBuffer.begin(), 
                                 energyBuffer.end(), 
                                 energy.begin());
                }
            }
            
            if(compToMeasure.force == true){
                uninitialized_cached_vector<real4> forceBuffer(this->pg->getNumberParticles());

                //Copy force to buffer
                {
                    auto force = this->pd->getForce(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 force.begin(), 
                                 force.end(), 
                                 forceBuffer.begin());
                }

                this->integrator->resetForce();
                this->integrator->updateForce();
                real3 totalForce = Measures::totalForce(pg, 
                                                        st); 
                
                //Copy back buffer to force
                {
                    auto force = this->pd->getForce(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 forceBuffer.begin(), 
                                 forceBuffer.end(), 
                                 force.begin());
                }
            }
            
            if(compToMeasure.virial == true){
                uninitialized_cached_vector<real> virialBuffer(this->pg->getNumberParticles());

                //Copy virial to buffer
                {
                    auto virial = this->pd->getVirial(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 virial.begin(), 
                                 virial.end(), 
                                 virialBuffer.begin());
                }

                this->integrator->resetVirial();
                this->integrator->updateVirial();
                real totalVirial = Measures::totalVirial(pg, 
                                                         st); 
                
                //Copy back buffer to virial
                {
                    auto virial = this->pd->getVirial(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 virialBuffer.begin(), 
                                 virialBuffer.end(), 
                                 virial.begin());
                }
            }
            
            if(compToMeasure.stress == true){
                uninitialized_cached_vector<tensor3> stressBuffer(this->pg->getNumberParticles());

                //Copy stress to buffer
                {
                    auto stress = this->pd->getStress(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 stress.begin(), 
                                 stress.end(), 
                                 stressBuffer.begin());
                }

                this->integrator->resetStress();
                this->integrator->updateStress();
                tensor3 totalStress = Measures::totalStress(pg, 
                                                            st); 
                
                //Copy back buffer to stress
                {
                    auto stress = this->pd->getStress(access::location::gpu, access::mode::read);     
                    thrust::copy(thrust::cuda::par.on(st),
                                 stressBuffer.begin(), 
                                 stressBuffer.end(), 
                                 stress.begin());
                }
            }
        }
};
    
template<class ForceField>
class EnergyMeasure: public SimulationStep{

        std::ofstream outPutFile;

        std::shared_ptr<ForceField> ff;

    public:
        
        EnergyMeasure(std::shared_ptr<ParticleGroup> pg,
                      int interval,
                      std::string outPutFileName,
                      std::shared_ptr<ForceField> ff):SimulationStep(pg,"EnergyMeasure",interval),
                                                      ff(ff){
            outPutFile = std::ofstream(outPutFileName);
        }

        void init(cudaStream_t st) override{}

        void applyStep(int step, cudaStream_t st) override{
            
            {
                auto energy = pd->getEnergy(access::location::gpu, access::mode::write);     
                thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), real(0));
            }

            uammd::Interactor::Computables comp;

            comp.energy = true;

            ff->sum(comp,st);

            real kE = Measures::totalKineticEnergy(this->pg,st);
            real pE = Measures::totalPotentialEnergy(this->pg,st);

            this->outPutFile << step << " " << std::setprecision(25) << kE << " " 
                                            << std::setprecision(25) << pE << " " 
                                            << std::setprecision(25) << kE+pE << std::endl;
        }
};

class InertiaMeasure: public SimulationStep{

        std::ofstream outPutFile;

    public:
        
        InertiaMeasure(std::shared_ptr<ParticleGroup> pg,
                       int interval,
                       std::string outPutFileName):SimulationStep(pg,"InertiaMeasure",interval){
            outPutFile = std::ofstream(outPutFileName);
        }

        void init(cudaStream_t st) override{}

        void applyStep(int step, cudaStream_t st) override{

            real totalMass = Measures::totalMass(this->pg,st);
            real3 com = Measures::centerOfMassPos(this->pg,totalMass,st);
            tensor3 I = Measures::inertiaTensor(this->pg,com,st);

            this->outPutFile << step << " " << I << std::endl;
        }
};

template<class Units_>
class TemperatureMeasure: public SimulationStep{

    private:

        std::ofstream outPutFile;
       
    public:

        TemperatureMeasure(std::shared_ptr<ParticleGroup> pg,
                           int interval,
                           std::string outPutFileName):SimulationStep(pg,"Temperature",interval){

            outPutFile = std::ofstream(outPutFileName);
        }
        
        void init(cudaStream_t st) override {}

        void applyStep(int step, cudaStream_t st) override{
            
            real kE = Measures::totalKineticEnergy(this->pg,st);
            
            int N = this->pg->getNumberParticles();

            real T = real(2.0/(3.0*N*Units_::KBOLTZ))*kE;

            this->outPutFile << step << " " << T << std::endl;
        }

};

class TemperaturePressureMeasure: public SimulationStep{

    private:

        real Kb;
        Box  box;

        std::ofstream outPutFile;
       
    public:

        TemperaturePressureMeasure(std::shared_ptr<ParticleGroup> pg,
                                   int interval,
                                   std::string outPutFileName,
                                   real Kb,Box box):SimulationStep(pg,"TemperaturePressure",interval),
                                                    Kb(Kb),
                                                    box(box){
            outPutFile = std::ofstream(outPutFileName);
        }
        
        void init(cudaStream_t st) override {}

        void applyStep(int step, cudaStream_t st) override{
            
            real kE = Measures::totalKineticEnergy(this->pg,st);
            
            tensor3 Pk = Measures::totalKineticPressure(this->pg,st);
            tensor3 Pv = Measures::totalStress(this->pg,st);

            int N = this->pg->getNumberParticles();

            real T = real(2.0/(3.0*N*Kb))*kE;

            real V = box.boxSize.x*
                     box.boxSize.y*
                     box.boxSize.z;

            real P = (Pk+Pv).trace()/(V*real(3.0));
            
            this->outPutFile << step << " " << T << " " << P << std::endl;
        }
        
        virtual void updateBox(Box newBox) override {
            box = newBox;
        }

};

template<class ForceField>
class StressMeasure: public SimulationStep{

        std::ofstream outPutFile;

        std::shared_ptr<ForceField> ff;

    public:
        
        StressMeasure(std::shared_ptr<ParticleGroup> pg,
                      int interval,
                      std::string outPutFileName,
                      std::shared_ptr<ForceField> ff):SimulationStep(pg,"StressMeasure",interval),
                                                      ff(ff){
            outPutFile = std::ofstream(outPutFileName);
        }

        void init(cudaStream_t st) override{}

        void applyStep(int step, cudaStream_t st) override{
            
            {
                auto stress = pd->getStress(access::location::gpu, access::mode::write);     
                thrust::fill(thrust::cuda::par.on(st), stress.begin(), stress.end(), real(0));
            }

            uammd::Interactor::Computables comp;

            comp.stress = true;

            ff->sum(comp,st);

            this->outPutFile << pg->getNumberParticles() << std::endl << step << std::endl;

            auto id   = pd->getId(access::location::cpu, 
                                  access::mode::read);

            auto pos   = pd->getPos(access::location::cpu, 
                                    access::mode::read);
            
            auto stress = pd->getStress(access::location::cpu, 
                                        access::mode::read);     

            auto groupIndex  = pg->getIndexIterator(access::location::cpu);
            auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            std::map<int,int> id_index;
            fori(0,pg->getNumberParticles()){
                int id_   = id[groupIndex[i]];
                int index = sortedIndex[id_];

                id_index[id_]=index;
            }

            for(const auto& ii : id_index){

                int index = ii.second;

                real4 pc = pos[index];
                real3 p = make_real3(pc);

                this->outPutFile << id[index]     << " " 
                                 << p             << " "
                                 << stress[index] << std::endl;
            }
        }
};

template<class NativeContactDefinition>
class NativeContactsMeasure: public SimulationStep{

    private:

        std::string fileNativeContacts;
        std::string labelNativeContacts;

        std::ofstream ncOutPutFile;
        std::ofstream comncOutPutFile;
        
        std::map<int,std::vector<int2>> molNativeContact;
        
        std::shared_ptr<NativeContactDefinition> ncd;

        std::map<int,std::shared_ptr<ParticleGroup>> mdlpg;
       
    public:

        NativeContactsMeasure(std::shared_ptr<ParticleGroup> pg,
                              int interval,
                              std::string fileNativeContacts,
                              std::string labelNativeContacts,
                              std::string ncOutPutFileName,
                              std::string comncOutPutFileName,
                              std::shared_ptr<NativeContactDefinition> ncd):SimulationStep(pg,"NaticeContactsMeasure",interval),
                                                                            fileNativeContacts(fileNativeContacts),
                                                                            labelNativeContacts(labelNativeContacts),
                                                                            ncd(ncd){

            ncOutPutFile    = std::ofstream(ncOutPutFileName);
            comncOutPutFile = std::ofstream(comncOutPutFileName);
        }
        
        void init(cudaStream_t st) override {
            
            std::shared_ptr<InputOutput::InputBlocksFile> nativeContactsDefinitionFile = std::make_shared<InputOutput::InputBlocksFile>(this->sys,
                                                                                                                                        fileNativeContacts);

            auto nativeContactBlock = nativeContactsDefinitionFile->getFileBlockIterator(labelNativeContacts);

            {
                auto molId  = pd->getModelId(access::location::cpu, access::mode::read);

                const int * sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                std::string line;
                std::stringstream parser;
                while (nativeContactBlock.next(line)){

                    parser.clear();
                    parser.str(line);

                    int id1;
                    int id2;

                    parser >> id1 >> id2;

                    if(parser.fail()){
                        sys->log<System::CRITICAL>("[NativeContacts] Error processing line \"%s\", while reading native contacts %s.", line.c_str(),labelNativeContacts.c_str());
                    }

                    int mol1 = molId[sortedIndex[id1]];
                    int mol2 = molId[sortedIndex[id2]];

                    if(mol1 == mol2){
                        molNativeContact[mol1].push_back({id1,id2});

                        sys->log<System::DEBUG1>("[NativeContacts] Added native contact (%i,%i), molecule (%i)",
                                id1,id2,
                                mol1);
                    } else {
                        sys->log<System::CRITICAL>("[NativeContacts] Error processing line \"%s\", particle ids (%i,%i) corresponds to different molecules (%i,%i)",
                                line.c_str(),id1,id2,
                                mol1,mol2);

                    }

                }

            }
            
            for(auto& molnc : molNativeContact){
                
                int m = molnc.first;
                selectors::modelId selector(m);
                std::shared_ptr<uammd::ParticleGroup> mpg = std::make_shared<uammd::ParticleGroup>(selector,pd,sys,"modelId_"+std::to_string(m));
                
                mdlpg[m]=mpg;

            }
        

        }

        void applyStep(int step, cudaStream_t st) override {


            std::map<int,int> molncCount;

            for(auto& molnc : molNativeContact){
                molncCount[molnc.first]=0;
            }
            
            for(auto& molnc : molNativeContact){
                for(int2& pair : molnc.second){
                    if(this->ncd->compute(pair.x,pair.y)){
                        molncCount[molnc.first]++;
                    }
                }
            }
            
            auto pos  = pd->getPos(access::location::cpu, access::mode::read);
            auto mass = pd->getMass(access::location::cpu, access::mode::read);
            
            comncOutPutFile << step << std::endl;
            ncOutPutFile << step << " ";
            for(auto const & ncc : molncCount){

                auto mpg = mdlpg[ncc.first];
                
                auto groupIndex  = mpg->getIndexIterator(access::location::cpu);

                real3 com = {0,0,0};
                
                real totalMass = 0;
                fori(0,mpg->getNumberParticles()){
                    int index = groupIndex[i];
                    totalMass+=mass[index];
                    com+=make_real3(pos[index])*mass[index];
                }

                comncOutPutFile << ncc.first  << " " << com/totalMass << " " << ncc.second << std::endl;
                ncOutPutFile    << ncc.second << " ";
            }
            ncOutPutFile << std::endl;;
        }

        virtual void updateBox(Box newBox) override {
            this->ncd->updateBox(newBox);
        }

};

}}

#endif
