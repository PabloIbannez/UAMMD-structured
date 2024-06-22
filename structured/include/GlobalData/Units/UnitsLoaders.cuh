#ifndef __UNITS_LOADER__
#define __UNITS_LOADER__
namespace uammd{
namespace structured{
namespace UnitsLoader{

    std::shared_ptr<typename Units::UnitsHandler>
    inline
    loadUnits(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string unitsType    = data.getType();
        std::string unitsSubType = data.getSubType();

        std::shared_ptr<typename Units::UnitsHandler> units;
        bool found = false;

        
        if("Units" == unitsType and "None" == unitsSubType){
            System::log<System::MESSAGE>("[UnitsLoader] (%s) Detected None units",path.back().c_str());
            units = std::make_shared<Units::None>(data);
            found = true;
        }
        if("Units" == unitsType and "KcalMol_A" == unitsSubType){
            System::log<System::MESSAGE>("[UnitsLoader] (%s) Detected KcalMol_A units",path.back().c_str());
            units = std::make_shared<Units::KcalMol_A>(data);
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[UnitsLoader] (%s) Could not find units %s::%s",path.back().c_str(),
                                           unitsType.c_str(),unitsSubType.c_str());
        }

        return units;

    }

    }}}
#endif
