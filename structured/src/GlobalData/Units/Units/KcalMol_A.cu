#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/UnitsFactory.cuh"

namespace uammd{
namespace structured{
namespace Units{

    class KcalMol_A: public UnitsHandler{

        public:

            KcalMol_A(DataEntry& data):UnitsHandler(data){}

            real getBoltzmannConstant()        override {return 1.987191e-03;}
            real getElectricConversionFactor() override {return 332.0716;}
    };

}}}

REGISTER_UNITS(
    Units,KcalMol_A,
    uammd::structured::Units::KcalMol_A
)
