#ifndef __UNITS_NONE__
#define __UNITS_NONE__

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

#endif
