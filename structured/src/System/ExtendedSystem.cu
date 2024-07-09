#include"System/ExtendedSystem.cuh"

namespace uammd{
namespace structured{

    void ExtendedSystem::loadSimulationInformation(std::string entryName){

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

    void ExtendedSystem::loadSimulationBackup(std::string entryName){

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

    void ExtendedSystem::init(){

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

    ExtendedSystem::ExtendedSystem(std::string inputFilePath):
    ExtendedSystem(0,nullptr,inputFilePath,{"system"}){}

    ExtendedSystem::ExtendedSystem(int argc, char *argv[],
                                   std::string inputFilePath,std::vector<std::string> path):
    uammd::System(argc,argv),
    path(path){
        //Check if the input file exists
        {
            std::ifstream file(inputFilePath);
            if(!file.good()){
                System::log<System::CRITICAL>("[ExtendedSystem] (%s) Input file not found: %s",
                                              path.back().c_str(),inputFilePath.c_str());
            }
        }

        try{
            input=std::make_shared<Input::Input>(inputFilePath);
        }catch(std::exception &e){
            System::log<System::CRITICAL>("[ExtendedSystem] (%s) Error reading input file: %s",
                                          path.back().c_str(),e.what());
        }

        ///////////////////////////////////////////

        this->init();

    }

    ///////////////////////////////////////////

    void ExtendedSystem::updateInputSystem(){

        System::log<System::DEBUG1>("[ExtendedSystem] (%s) Updating input backup",this->path.back().c_str());

        auto data = this->input->getDataEntry(backupPath);

        data.template setParameter<ullint>("seed",seed);

        if(saveBackup){
            data.template setParameter<bool>("restartedFromBackup",restartedFromBackup);
            data.template setParameter<ullint>("lastBackupStep",lastBackupStep);
        }

    }

    void ExtendedSystem::finish(){

        std::string line;
        fori(0,29) line += "━ ";
        System::log<System::MESSAGE>("%s", line.c_str());

        if(cudaStreamCreated){
            System::log<System::MESSAGE>("[ExtendedSystem] (%s) Destroying cuda stream",path.back().c_str());
            cudaStreamDestroy(stream);
        }
        uammd::System::finish();
    }

    std::shared_ptr<ExtendedSystem> getExtendedSystem(std::shared_ptr<uammd::System> sys){
        return std::static_pointer_cast<ExtendedSystem>(sys);
    }
}}
