#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/Units/KcalMol_A.cuh"

namespace uammd{
namespace structured{
namespace Units{

    KcalMol_A::KcalMol_A(DataEntry& data): UnitsHandler(data){}

    real KcalMol_A::getBoltzmannConstant(){
        return 1.987191E-03; // kcal/(mol*K)
    }

    real KcalMol_A::getElectricConversionFactor(){
        return 332.0716; // kcal/(mol*A)
    }

}}}
