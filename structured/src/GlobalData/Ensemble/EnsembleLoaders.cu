#include "GlobalData/Ensemble/EnsembleLoaders.cuh"

namespace uammd{
namespace structured{
namespace EnsembleLoader{

    bool isEnsembleAvailable(std::shared_ptr<ExtendedSystem> sys,
                             std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string ensembleType    = data.getType();
        std::string ensembleSubType = data.getSubType();

        return Ensemble::EnsembleFactory::getInstance().isEnsembleRegistered(ensembleType,ensembleSubType);
    }


    std::shared_ptr<typename Ensemble::EnsembleHandler>
    loadEnsemble(std::shared_ptr<ExtendedSystem> sys,
                 std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string ensembleType    = data.getType();
        std::string ensembleSubType = data.getSubType();

        //////////////////////////////////////////////////////////////

        return Ensemble::EnsembleFactory::getInstance().createEnsemble(ensembleType,ensembleSubType,
                                                                       data);
    }

}}}
