#ifndef __SIMULATION_STEP_OUTPUT__
#define __SIMULATION_STEP_OUTPUT__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationOutput{

class WriteStep: public SimulationStepBase{

    public:

  std::vector<std::string> outputFormatAvail={"coord",
                                              "sp","spo",
                                              "spf",
                                              "xyz",
                                              "pdb",
                                              "itpv","itpd",
                                              "dcd",
                                              "lammpstrj",
                                              "vel",
                                              "magnet", "xyzm", "spm",
					      "svv", "svvm", "svvma"};

    private:

        std::string outputFilePath;
        std::string outputFormat;

        bool pbc;
        bool append;

        /////////////////////////

        std::ofstream outputFile;

        int frame;

        Box getBox(){
            Box box = gd->getEnsemble()->getBox();
            box.setPeriodicity(pbc,pbc,pbc);
            return box;
        }

    public:

        WriteStep(std::shared_ptr<ParticleGroup>  pg,
                  std::shared_ptr<IntegratorManager> integrator,
                  std::shared_ptr<ForceField>    ff,
                  DataEntry& data,
                  std::string name):SimulationStepBase(pg,integrator,ff,data,name){

            //Read parameters

            outputFilePath = data.getParameter<std::string>("outputFilePath");
            outputFormat   = data.getParameter<std::string>("outputFormat");

            pbc    = data.getParameter<bool>("pbc",true);
            append = data.getParameter<bool>("append",false);

            //Print parameters
            System::log<System::MESSAGE>("[WriteStep] (%s) outputFilePath: %s",name.c_str(),outputFilePath.c_str());
            System::log<System::MESSAGE>("[WriteStep] (%s) outputFormat: %s",name.c_str(),outputFormat.c_str());

            System::log<System::MESSAGE>("[WriteStep] (%s) pbc: %i",name.c_str(),pbc);
            System::log<System::MESSAGE>("[WriteStep] (%s) append: %i",name.c_str(),append);

        }

        void init(cudaStream_t st) override{

            //Add backup

            if(true){
                if(this->intervalStep != 0){
                    if(std::find(outputFormatAvail.begin(),
                                 outputFormatAvail.end(),
                                 outputFormat)==outputFormatAvail.end()){

                        std::string message = "[WriteStep] (" + this->name + ") The selected output format ("+outputFormat+") is not available,options are:";

                        int i;
                        for(i=0;i<outputFormatAvail.size()-1;i++){
                            message.append(" "+outputFormatAvail[i]+",");
                        }
                        message.append(" "+outputFormatAvail[i]+".");

                        System::log<System::CRITICAL>("%s", message.c_str());
                    } else {

                        if(outputFormat ==
                           std::string("dcd")){

                            if(!append){
                                std::ofstream outPSF(outputFilePath + ".psf");

                                psf::WritePSF(pg,outPSF);

                                outputFilePath = outputFilePath + "." + outputFormat;
                                outputFile.open(outputFilePath,std::ios::binary);

                                InputOutput::Output::WriteDCDheader(pg,
                                                                    frame,
                                                                    this->intervalStep,
                                                                    outputFile);
                            } else {
                                outputFilePath = outputFilePath + "." + outputFormat;
                                outputFile.open(outputFilePath,std::ios::binary | std::ios_base::app);
                            }
                        } else {
                            outputFilePath = outputFilePath + "." + outputFormat;
                            if(!append){
                                outputFile.open(outputFilePath);
                            } else {
                                outputFile.open(outputFilePath,std::ios_base::app);
                            }
                        }
                    }
                } else {
                    System::log<System::WARNING>("[WriteStep] (%s) Write interval is 0, no files will be generated!",name.c_str());
                }
            }
        }

        void applyStep(ullint step, cudaStream_t st) override{

            if(step==this->getLastStepApplied()){
                return;
            }

            frame++;

            if(outputFormat ==
               std::string("coord")){
                InputOutput::Output::WriteCoord(pg,this->getBox(),outputFile);
            } else if(outputFormat ==
                      std::string("sp")){
                InputOutput::Output::WriteSP(pg,this->getBox(),outputFile);
            } else if(outputFormat ==
                      std::string("spo")){
                InputOutput::Output::WriteSPO(pg,this->getBox(),outputFile);
            } else if(outputFormat ==
                      std::string("spf")){
                InputOutput::Output::WriteSPF(pg,this->getBox(),outputFile);
            } else if (outputFormat ==
                       std::string("xyz")){
                InputOutput::Output::WriteXYZ(pg,this->getBox(),outputFile);
            } else if (outputFormat ==
                       std::string("pdb")){
                InputOutput::Output::WritePDB(pg,this->getBox(),frame,outputFile);
            } else if (outputFormat ==
                       std::string("itpv")){
                InputOutput::Output::WriteITPV(pg,this->getBox(),outputFile);
            } else if (outputFormat ==
                       std::string("itpd")){
                InputOutput::Output::WriteITPD(pg,this->getBox(),outputFile);
            } else if (outputFormat ==
                       std::string("dcd")){
                InputOutput::Output::WriteDCD(pg,
                                              this->getBox(),
                                              frame,step,
                                              outputFile);
            } else if (outputFormat ==
                       std::string("lammpstrj")){
                real time = step*gd->getFundamental()->getTimeStep();
                InputOutput::Output::WriteLAMMPS(pg,this->getBox(),time,outputFile);
            } else if (outputFormat ==
                       std::string("vel")){
                InputOutput::Output::WriteVelocity(pg,outputFile);
	    } else if (outputFormat ==
                       std::string("magnet")){
	      InputOutput::Output::WriteMagnetization(pg,outputFile);
            } else if (outputFormat ==
                       std::string("xyzm")){
	      InputOutput::Output::WriteXYZMagnetization(pg,outputFile);
            } else if (outputFormat ==
                       std::string("spm")){
	      InputOutput::Output::WriteSPM(pg,this->getBox(),outputFile);
	    } else if (outputFormat ==
                       std::string("svv")){
	      InputOutput::Output::WriteSVV(pg,this->getBox(),outputFile);
	    } else if (outputFormat ==
                       std::string("svvm")){
	      InputOutput::Output::WriteSVVM(pg,this->getBox(),outputFile);
	    } else if (outputFormat ==
                       std::string("svvma")){
	      InputOutput::Output::WriteSVVMA(pg,this->getBox(),outputFile);
            } else {
	      System::log<System::CRITICAL>("[WriteStep] (%s) It should not have happen,"
					    "have the output format changed ?",name.c_str());
            }
        }
};

class WritePatchyParticlesStep: public WriteStep{

        std::shared_ptr<GlobalData>   gd;

        std::shared_ptr<ParticleData>  pdParent;
        std::shared_ptr<ParticleGroup> pgParent;

        std::map<std::string,std::shared_ptr<Interactor::PatchyParticles_<>>> patches;
        std::map<std::string,std::shared_ptr<ParticleGroup>> patchesGroups;

        std::shared_ptr<ParticleData> pdBuffer;

        void updateParticleDataBuffer(cudaStream_t st){

            for(auto& p : patches){
                p.second->updatePatchyParticles(st);
            }

            //Update the particle data buffer
            auto posBuffer    = pdBuffer->getPos(access::location::cpu,access::mode::write);
            auto radiusBuffer = pdBuffer->getRadius(access::location::cpu,access::mode::write);

            //Load parent particle data
            int offset = 0;
            {
                auto pos    = pdParent->getPos(access::location::cpu,access::mode::read);
                auto radius = pdParent->getRadius(access::location::cpu,access::mode::read);

                auto groupIndex = this->pgParent->getIndexIterator(access::location::cpu);

                for(int i=0;i<this->pgParent->getNumberParticles();i++){
                    int index = groupIndex[i];
                    posBuffer[offset]    = pos[index];
                    radiusBuffer[offset] = radius[index];
                    offset++;
                }
            }
            int typesOffset = gd->getTypes()->getNumberOfTypes()+1;

            //Load patch particle data
            for(auto& p : patches){

                auto pos    = p.second->getPatchesParticleData()->getPos(access::location::cpu,access::mode::read);
                auto radius = p.second->getPatchesParticleData()->getRadius(access::location::cpu,access::mode::read);

                auto groupIndex = patchesGroups[p.first]->getIndexIterator(access::location::cpu);

                for(int i=0;i<patchesGroups[p.first]->getNumberParticles();i++){
                    int index = groupIndex[i];
                    posBuffer[offset]    = pos[index];
                    posBuffer[offset].w  = int(posBuffer[offset].w) + typesOffset;
                    radiusBuffer[offset] = radius[index];
                    offset++;
                }

                typesOffset += p.second->getPatchesGlobalData()->getTypes()->getNumberOfTypes()+1;
            }
        };

    public:

        WritePatchyParticlesStep(std::shared_ptr<ParticleGroup>     pg,
                                 std::shared_ptr<IntegratorManager> integrator,
                                 std::shared_ptr<ForceField>        ff,
                                 DataEntry& data,
                                 std::string name):WriteStep(pg,integrator,ff,data,name){

            //Store the global data
            gd = topology->getGlobalData();
            //Store the original particle data
            pdParent = pg->getParticleData();
            pgParent = pg;

            std::map<std::string,std::shared_ptr<typename uammd::Interactor>> patchesTmp = topology->getInteractorsByClass("PatchyParticles");

            for(auto& p : patchesTmp){
                patches[p.first]=std::static_pointer_cast<Interactor::PatchyParticles_<>>(p.second);
                System::log<System::MESSAGE>("[WritePatchyParticlesStep] (%s) Found PatchyParticles interactor \"%s\"",name.c_str(),p.first.c_str());
            }

            if(patches.size()==0){
                System::log<System::WARNING>("[WritePatchyParticlesStep] (%s) No PatchyParticles interactor found",name.c_str());
            }
        }

        void init(cudaStream_t st) {

            //Get the ids of the particles in the group
            auto pgIds = pgParent->getPropertyIterator(pdParent->getId(access::location::cpu, access::mode::read).begin(),
                                                       access::location::cpu);

            //Iterate over all PatchyParticles interactors and compute the number of particles
            int Np = pgParent->getNumberParticles();
            for(auto p : patches){
                auto patchesId = p.second->getPatchesParticleData()->getId(access::location::cpu,access::mode::read);
                auto parentId  = p.second->getPatchesParticleData()->getModelId(access::location::cpu,access::mode::read);

                //Find the id of the patches whose parent is in pgIds
                std::vector<int> ids;
                for(int i=0;i<pgParent->getNumberParticles();i++){
                    for(int j=0;j<p.second->getPatchesParticleData()->getNumParticles();j++){
                        if(pgIds[i]==parentId[j]){
                            ids.push_back(patchesId[j]);
                        }
                    }
                }

                //Create a new group with the patches
                GroupUtils::selectors::ids selector(ids);
                patchesGroups[p.first] = std::make_shared<ParticleGroup>(selector,
                                                                         p.second->getPatchesParticleData(),
                                                                         "WritePatchyParticlesStep_"+p.first+pgParent->getName());

                Np += ids.size();
            }

            pdBuffer = std::make_shared<ParticleData>(Np,pgParent->getParticleData()->getSystem());
            //Overriding the particle group
            this->pg = std::make_shared<ParticleGroup>(pdBuffer,"ParentPatches_"+pgParent->getName());

            WriteStep::init(st);
        }

        void applyStep(ullint step, cudaStream_t st){
            this->updateParticleDataBuffer(st);
            WriteStep::applyStep(step,st);
        }
};

}}}}

#endif
