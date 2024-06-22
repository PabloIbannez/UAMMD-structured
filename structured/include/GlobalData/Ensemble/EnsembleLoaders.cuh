#ifndef __ENSEMBLE_LOADER__
#define __ENSEMBLE_LOADER__
namespace uammd{
namespace structured{
namespace EnsembleLoader{

    std::shared_ptr<typename Ensemble::EnsembleHandler>
    inline
    loadEnsemble(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string ensembleType    = data.getType();
        std::string ensembleSubType = data.getSubType();

        std::shared_ptr<typename Ensemble::EnsembleHandler> ensemble;
        bool found = false;

        
        if("Ensemble" == ensembleType and "NVT" == ensembleSubType){
            System::log<System::MESSAGE>("[EnsembleLoader] (%s) Detected NVT ensemble",path.back().c_str());
            ensemble = std::make_shared<Ensemble::NVT>(data);
            found = true;
        }
        if("Ensemble" == ensembleType and "NVTlambda" == ensembleSubType){
            System::log<System::MESSAGE>("[EnsembleLoader] (%s) Detected NVTlambda ensemble",path.back().c_str());
            ensemble = std::make_shared<Ensemble::NVTlambda>(data);
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[EnsembleLoader] (%s) Could not find ensemble %s::%s",path.back().c_str(),
                                           ensembleType.c_str(),ensembleSubType.c_str());
        }

        return ensemble;

    }

    }}}
#endif
