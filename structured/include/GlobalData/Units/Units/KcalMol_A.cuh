#ifndef __UNITS_KCALMOL_A__
#define __UNITS_KCALMOL_A__

namespace uammd{
namespace structured{
namespace Units{

    class KcalMol_A: public UnitsHandler{

        public:

            KcalMol_A(DataEntry& data);

            real getBoltzmannConstant()        override;
            real getElectricConversionFactor() override;
    };

}}}

#endif
