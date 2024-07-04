#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/UnitsFactory.cuh"

namespace uammd{
namespace structured{
namespace Units{

    class None: public UnitsHandler{

        public:

            None(DataEntry& data):UnitsHandler(data){}

            real getBoltzmannConstant()        override {return 1.0;}
            real getElectricConversionFactor() override {return 1.0;}
    };

}}}

REGISTER_UNITS(
    Units,None,
    uammd::structured::Units::None
)
