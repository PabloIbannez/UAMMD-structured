#ifndef __TYPES_LOADER__
#define __TYPES_LOADER__
namespace uammd{
namespace structured{
namespace TypesLoader{

    std::shared_ptr<typename Types::TypesHandler>
    inline
    loadTypes(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string typesType    = data.getType();
        std::string typesSubType = data.getSubType();

        std::shared_ptr<typename Types::TypesHandler> types;
        bool found = false;

        
        if("Types" == typesType and "Basic" == typesSubType){
            System::log<System::MESSAGE>("[TypesLoader] (%s) Detected Basic types",path.back().c_str());
            types = std::make_shared<Types::Basic>(data);
            found = true;
        }
        if("Types" == typesType and "Polarizable" == typesSubType){
            System::log<System::MESSAGE>("[TypesLoader] (%s) Detected Polarizable types",path.back().c_str());
            types = std::make_shared<Types::Polarizable>(data);
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[TypesLoader] (%s) Could not find types %s::%s",path.back().c_str(),
                                           typesType.c_str(),typesSubType.c_str());
        }

        return types;

    }

    }}}
#endif
