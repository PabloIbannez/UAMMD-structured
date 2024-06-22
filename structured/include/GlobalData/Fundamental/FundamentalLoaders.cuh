#ifndef __FUNDAMENTAL_LOADER__
#define __FUNDAMENTAL_LOADER__
namespace uammd{
namespace structured{
namespace FundamentalLoader{

    std::shared_ptr<typename Fundamental::FundamentalHandler>
    inline
    loadFundamental(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string fundamentalType    = data.getType();
        std::string fundamentalSubType = data.getSubType();

        std::shared_ptr<typename Fundamental::FundamentalHandler> fundamental;
        bool found = false;

        
        if("Fundamental" == fundamentalType and "Time" == fundamentalSubType){
            System::log<System::MESSAGE>("[FundamentalLoader] (%s) Detected Time fundamental",path.back().c_str());
            fundamental = std::make_shared<Fundamental::Time>(data);
            found = true;
        }
        if("Fundamental" == fundamentalType and "DynamicallyBondedPatchyParticles" == fundamentalSubType){
            System::log<System::MESSAGE>("[FundamentalLoader] (%s) Detected DynamicallyBondedPatchyParticles fundamental",path.back().c_str());
            fundamental = std::make_shared<Fundamental::DynamicallyBondedPatchyParticles>(data);
            found = true;
        }
        if("Fundamental" == fundamentalType and "None" == fundamentalSubType){
            System::log<System::MESSAGE>("[FundamentalLoader] (%s) Detected None fundamental",path.back().c_str());
            fundamental = std::make_shared<Fundamental::None>(data);
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[FundamentalLoader] (%s) Could not find fundamental %s::%s",path.back().c_str(),
                                           fundamentalType.c_str(),fundamentalSubType.c_str());
        }

        return fundamental;

    }

    }}}
#endif
