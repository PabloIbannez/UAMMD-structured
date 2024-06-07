#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/Units/None.cuh"

namespace uammd{
namespace structured{
namespace Units{

    None::None(DataEntry& data): UnitsHandler(data){}

    real None::getBoltzmannConstant(){
        return 1.0;
    }

    real None::getElectricConversionFactor(){
        return 1.0;
    }

}}}
