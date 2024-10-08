#include "GlobalData/GlobalData.cuh"

namespace uammd{
namespace structured{

    GlobalDataBase::GlobalDataBase(std::shared_ptr<ExtendedSystem> sys,
            std::vector<std::string> path):sys(sys),path(path){

        std::shared_ptr<InputEntryManager>   globalInfo;
        globalInfo = std::make_shared<InputEntryManager>(sys,path);

        //Units
        std::vector<std::string> unitsPath = path;
        unitsPath.push_back("units");
        if(globalInfo->isEntryPresent("units")){
            unitsHandler = UnitsLoader::loadUnits(sys,unitsPath);
        }else{
            System::log<System::MESSAGE>("[GlobalDataBase] Units not specified, using default units, \"None\"");
            auto input = sys->getInput();
            input->addEntry(unitsPath,"Units","None");

            unitsHandler = UnitsLoader::loadUnits(sys,unitsPath);
        }

        //Ensemble
        std::vector<std::string> ensemblePath = path;
        ensemblePath.push_back("ensemble");

        ensembleHandler = EnsembleLoader::loadEnsemble(sys,ensemblePath);

        System::log<System::DEBUG1>("[GlobalDataBase] GlobalDataBase created");
    }

    void GlobalDataBase::updateInputGlobalData(){

        //Ensemble
        std::vector<std::string> ensemblePath = path;
        ensemblePath.push_back("ensemble");

        ensembleHandler->updateDataEntry(sys->getInput()->getDataEntry(ensemblePath));
    }

    //GlobalData

    void GlobalData::init(){

        std::shared_ptr<ExtendedSystem> sys = globalDataBase->getSystem();

        std::shared_ptr<InputEntryManager>   globalInfo;
        globalInfo = std::make_shared<InputEntryManager>(sys,path);

        //Fundamental
        std::vector<std::string> fundamentalPath = path;
        fundamentalPath.push_back("fundamental");
        if(globalInfo->isEntryPresent("fundamental")){
            fundamentalHandler = FundamentalLoader::loadFundamental(sys,fundamentalPath);
        }else{
            System::log<System::MESSAGE>("[GlobalDataBase] Fundamental not specified, using default fundamental, \"Time\"");
            auto input = sys->getInput();
            input->addEntry(fundamentalPath,"Fundamental","Time");

            fundamentalHandler = FundamentalLoader::loadFundamental(sys,fundamentalPath);
        }

        //Load types
        std::vector<std::string> typesPath = path;
        typesPath.push_back("types");

        typesHandler = TypesLoader::loadTypes(sys,typesPath);

    }


    GlobalData::GlobalData(std::shared_ptr<GlobalDataBase>  globalDataBase,
                           std::vector<std::string>         path):globalDataBase(globalDataBase),path(path){
        this->init();
    }

    GlobalData::GlobalData(std::shared_ptr<ExtendedSystem>  sys,
                           std::vector<std::string>         path):path(path){

        globalDataBase = std::make_shared<GlobalDataBase>(sys,path);

        this->init();

    }

    GlobalData::GlobalData(std::shared_ptr<ExtendedSystem>  sys):GlobalData(sys,{"global"}){}

    void GlobalData::updateInputGlobalData(){

        std::shared_ptr<ExtendedSystem> sys = globalDataBase->getSystem();

        //Fundamental
        std::vector<std::string> fundamentalPath = path;
        fundamentalPath.push_back("fundamental");

        fundamentalHandler->updateDataEntry(sys->getInput()->getDataEntry(fundamentalPath));

        globalDataBase->updateInputGlobalData();
    }

}}
