#ifndef __UNITS_NONE__
#define __UNITS_NONE__

namespace uammd{
namespace structured{
namespace Units{

    class None: public UnitsHandler{

        public:

            None(DataEntry& data);

            real getBoltzmannConstant()        override;
            real getElectricConversionFactor() override;
    };

}}}

#endif
