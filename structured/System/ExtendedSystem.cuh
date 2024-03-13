#ifndef __EXTENDED_SYSTEM__
#define __EXTENDED_SYSTEM__

#include"ThirdParty/json.hpp"

namespace uammd{
namespace structured{

    template<class InputType_>
    class ExtendedSystem_ : public uammd::System {

        public:

            using InputType = InputType_;

            //Enum available states, running and stopped
            enum SIMULATION_STATE {RUNNING, STOPPED};

        private:

            //Simulation state
            SIMULATION_STATE state = RUNNING;

            bool cudaStreamCreated = false;
            cudaStream_t stream;

            std::shared_ptr<InputType> input;

            ///////////////////////////
            std::vector<std::string> path;

            std::vector<std::string> infoPath;
            std::vector<std::string> backupPath;

            std::string name;
            ullint seed;

            //Backup
            bool        saveBackup;
            std::string backupFilePath;

            bool    restartedFromBackup;
            ullint  backupRestartStep;

            ullint lastBackupStep;

            ullint backupStartStep;
            ullint backupEndStep;
            ullint backupIntervalStep;

            void loadSimulationInformation(std::string entryName){

                std::vector<std::string> infoPath = path;
                infoPath.push_back(entryName);

                auto data = input->getDataEntry(infoPath);

                name = data.template getParameter<std::string>("name");

                if(data.isParameterAdded("seed")){
                        seed = data.template getParameter<ullint>("seed");
                } else {
                        seed = std::chrono::system_clock::now().time_since_epoch().count();
                }

                System::log<System::MESSAGE>("[ExtendedSystem] (%s) Name: %s",this->path.back().c_str(),name.c_str());
                System::log<System::MESSAGE>("[ExtendedSystem] (%s) Seed: %llu",this->path.back().c_str(),seed);


            }

            void loadSimulationBackup(std::string entryName){

                backupPath = path;
                backupPath.push_back(entryName);

                auto data = input->getDataEntry(backupPath);

                //Backup
                restartedFromBackup = data.template getParameter<bool>("restartedFromBackup",false);

                backupFilePath = data.template getParameter<std::string>("backupFilePath");

                lastBackupStep     = data.template getParameter<ullint>("lastBackupStep",0);

                backupStartStep    = data.template getParameter<ullint>("backupStartStep",0);
                backupEndStep      = data.template getParameter<ullint>("backupEndStep",std::numeric_limits<ullint>::max());
                backupIntervalStep = data.template getParameter<ullint>("backupIntervalStep");

                if(restartedFromBackup){
                    backupRestartStep = lastBackupStep;
                }
            }

            void init(){

                this->saveBackup = false;
                this->restartedFromBackup = false;

                bool informationLoaded = false;
                for(std::string entryName : input->getEntriesList(this->path)){
                    std::vector<std::string> entryPath = this->path;
                    entryPath.push_back(entryName);

                    auto dataEntry = input->getDataEntry(entryPath);

                    std::string entryType    = dataEntry.getType();
                    std::string entrySubType = dataEntry.getSubType();

                    if(entryType == "Simulation"){
                        if(entrySubType == "Information"){
                            loadSimulationInformation(entryName);
                            informationLoaded = true;
                        }
                        else if(entrySubType == "Backup"){
                            saveBackup = true;
                            loadSimulationBackup(entryName);
                        } else {
                            System::log<System::CRITICAL>("[ExtendedSystem] (%s) Unknown Simulation entry type, subType: %s, %s",
                                                           this->path.back().c_str(),entryType.c_str(),entrySubType.c_str());
                        }
                    } else {
                        System::log<System::CRITICAL>("[ExtendedSystem] (%s) Unknown entry type: %s",
                                                      this->path.back().c_str(),
                                                      entryType.c_str());

                    }
                }

                if(!informationLoaded){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) No Simulation Information found!",
                                                   this->path.back().c_str());
                }

            }

        public:

            ExtendedSystem_(std::string inputFilePath):ExtendedSystem_(0,nullptr,inputFilePath,{"system"}){}
            ExtendedSystem_(int argc, char *argv[],std::string inputFilePath,std::vector<std::string> path):
                                                                              uammd::System(argc,argv),
                                                                              path(path){
                //Check if the input file exists
                {
                    std::ifstream file(inputFilePath);
                    if(!file.good()){
                        System::log<System::CRITICAL>("[ExtendedSystem] (%s) Input file not found: %s",path.back().c_str(),inputFilePath.c_str());
                    }
                }

                try{
                    input=std::make_shared<InputType>(inputFilePath);
                }catch(std::exception &e){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Error reading input file: %s",path.back().c_str(),e.what());
                }

                ///////////////////////////////////////////

                this->init();

            }

            cudaStream_t getCudaStream(){
                if(!cudaStreamCreated){
                    System::log<System::MESSAGE>("[ExtendedSystem] (%s) Creating cuda stream",path.back().c_str());
                    CudaSafeCall(cudaStreamCreate(&stream));
                    cudaStreamCreated = true;
                }
                return stream;
            }

            //void synchronizeCudaStream(bool force=false){

            //    if(!cudaStreamCreated){
            //        System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to synchronize a non created stream",path.back().c_str());
            //    }

            //    if(force){
            //        cudaStreamSynchronize(stream);
            //    } else {
            //        cudaStreamSynchronize(stream);
            //    }
            //}

            SIMULATION_STATE getState(){
                return state;
            }

            void setState(SIMULATION_STATE state){
                this->state = state;
            }

            std::shared_ptr<InputType> getInput(){
                return input;
            }

            //Getters
            std::string getName(){ return name; }
            ullint getSeed(){ return seed; }

            //Setters
            void setSeed(ullint seed){ this->seed = seed; }

            //Backup getters
            bool getSaveBackup(){ return saveBackup; }

            std::string getBackupFilePath(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup file path,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupFilePath;
            }

            bool getRestartedFromBackup(){ return restartedFromBackup; }

            ullint getBackupRestartStep(){
                if(!restartedFromBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup restart step, "
                                                  "but the simulation is not a backup!",path.back().c_str());
                }
                return backupRestartStep;
            }

            ullint getLastBackupStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get last backup step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return lastBackupStep;
            }

            ullint getBackupStartStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup start step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupStartStep;
            }

            ullint getBackupEndStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup end step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupEndStep;
            }

            ullint getBackupIntervalStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup interval step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupIntervalStep;
            }

            //Backup setters

            void setRestartedFromBackup(bool restartedFromBackup){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set restartedFromBackup,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                this->restartedFromBackup = restartedFromBackup;
            }

            void setBackupRestartStep(ullint backupRestartStep){
                if(!restartedFromBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set backup restart step,"
                                                  " but the simulation is not a backup!",path.back().c_str());
                }
                this->backupRestartStep = backupRestartStep;
            }

            void setLastBackupStep(ullint lastBackupStep){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set last backup step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                this->lastBackupStep = lastBackupStep;
            }

            void updateInputSystem(){

                System::log<System::DEBUG1>("[ExtendedSystem] (%s) Updating input backup",this->path.back().c_str());

                auto data = this->input->getDataEntry(backupPath);

                data.template setParameter<ullint>("seed",seed);

                if(saveBackup){
                    data.template setParameter<bool>("restartedFromBackup",restartedFromBackup);
                    data.template setParameter<ullint>("lastBackupStep",lastBackupStep);
                }

            }

            ///////////////////////////////////////////

            void finish(){

                std::string line;
                fori(0,29) line += "━ ";
                System::log<System::MESSAGE>("%s", line.c_str());

                if(cudaStreamCreated){
                    System::log<System::MESSAGE>("[ExtendedSystem] (%s) Destroying cuda stream",path.back().c_str());
                    cudaStreamDestroy(stream);
                }
                uammd::System::finish();
            }
    };

    using ExtendedSystem = ExtendedSystem_<InputJSON::InputJSON>;

    using ParametersEntry = ExtendedSystem::InputType::parametersEntry;
    using TypeEntry       = ExtendedSystem::InputType::typeEntry;
    using DataEntry       = ExtendedSystem::InputType::dataEntry;

    std::shared_ptr<ExtendedSystem> getExtendedSystem(std::shared_ptr<uammd::System> sys){
        return std::static_pointer_cast<ExtendedSystem>(sys);
    }
}}

#endif
