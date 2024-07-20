#include "GlobalData/Types/TypesLoaders.cuh"

namespace uammd{
namespace structured{
namespace TypesLoader{

    bool isTypesAvailable(std::shared_ptr<ExtendedSystem> sys,
                          std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string typesType    = data.getType();
        std::string typesSubType = data.getSubType();

        return Types::TypesFactory::getInstance().isTypesRegistered(typesType,typesSubType);
    }


    std::shared_ptr<typename Types::TypesHandler>
    loadTypes(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string typesType    = data.getType();
        std::string typesSubType = data.getSubType();

        //////////////////////////////////////////////////////////////

        return Types::TypesFactory::getInstance().createTypes(typesType,typesSubType,
                                                              data);
    }

}}}
