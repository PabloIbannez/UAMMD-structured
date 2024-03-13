#ifndef __UNITS_KCALMOL_A__
#define __UNITS_KCALMOL_A__

namespace uammd{
namespace structured{
namespace Units{

    class KcalMol_A: public UnitsHandler{

        public:

            KcalMol_A(DataEntry& data):UnitsHandler(data){}

            real getBoltzmannConstant()        override {return 1.987191E-03;}
            real getElectricConversionFactor() override {return 332.0716;}
    };

}}}

#endif
