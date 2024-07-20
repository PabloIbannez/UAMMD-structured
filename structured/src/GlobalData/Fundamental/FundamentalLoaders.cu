#include "GlobalData/Fundamental/FundamentalLoaders.cuh"

namespace uammd{
namespace structured{
namespace FundamentalLoader{

    bool isFundamentalAvailable(std::shared_ptr<ExtendedSystem> sys,
                                std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string fundamentalType    = data.getType();
        std::string fundamentalSubType = data.getSubType();

        return Fundamental::FundamentalFactory::getInstance().isFundamentalRegistered(fundamentalType,fundamentalSubType);
    }


    std::shared_ptr<typename Fundamental::FundamentalHandler>
    loadFundamental(std::shared_ptr<ExtendedSystem> sys,
                    std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string fundamentalType    = data.getType();
        std::string fundamentalSubType = data.getSubType();

        //////////////////////////////////////////////////////////////

        return Fundamental::FundamentalFactory::getInstance().createFundamental(fundamentalType,fundamentalSubType,
                                                                                data);
    }

}}}
