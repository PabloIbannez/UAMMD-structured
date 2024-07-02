#include"System/ExtendedSystem.cuh"

namespace uammd{
namespace structured{

    template<class InputType_>
    void ExtendedSystem_<InputType_>::loadSimulationInformation(std::string entryName){

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

    template<class InputType_>
    void ExtendedSystem_<InputType_>::loadSimulationBackup(std::string entryName){

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

    template<class InputType_>
    void ExtendedSystem_<InputType_>::init(){

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

    template<class InputType_>
    ExtendedSystem_<InputType_>::ExtendedSystem_(std::string inputFilePath):
    ExtendedSystem_<InputType_>(0,nullptr,inputFilePath,{"system"}){}

    template<class InputType_>
    ExtendedSystem_<InputType_>::ExtendedSystem_(int argc, char *argv[],
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
            input=std::make_shared<InputType>(inputFilePath);
        }catch(std::exception &e){
            System::log<System::CRITICAL>("[ExtendedSystem] (%s) Error reading input file: %s",
                                          path.back().c_str(),e.what());
        }

        ///////////////////////////////////////////

        this->init();

    }

    ///////////////////////////////////////////

    template<class InputType_>
    void ExtendedSystem_<InputType_>::updateInputSystem(){

        System::log<System::DEBUG1>("[ExtendedSystem] (%s) Updating input backup",this->path.back().c_str());

        auto data = this->input->getDataEntry(backupPath);

        data.template setParameter<ullint>("seed",seed);

        if(saveBackup){
            data.template setParameter<bool>("restartedFromBackup",restartedFromBackup);
            data.template setParameter<ullint>("lastBackupStep",lastBackupStep);
        }

    }

    template<class InputType_>
    void ExtendedSystem_<InputType_>::finish(){

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

    template class ExtendedSystem_<InputJSON::InputJSON>;
}}
