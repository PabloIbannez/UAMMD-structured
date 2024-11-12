#include "GlobalData/Units/UnitsHandler.cuh"
#include "GlobalData/Units/UnitsFactory.cuh"

namespace uammd{
namespace structured{
namespace Units{

    class lenSI: public UnitsHandler{
        /*
           International System units, have a parameter lenghtScale to scale the Lenght.

           length      = lenghtScale  m
           mass        = lenghScale^3 g  (because mass scale with volume)
           time        = lenghtScale  s  (So energy scale with volume)
           charge      = lenghtScale  C  (So voltage is in Volts)
           Temperature = Kelvin

           */
            real scale3;
            real scale2;


        public:

            lenSI(DataEntry& data):UnitsHandler(data){

                real scale  = data.getParameter<real>("lengthScale",1.0);
                scale3 = scale*scale*scale;
                scale2 = scale*scale;
            }

            real getBoltzmannConstant()        override {return 1.380649e-23/scale3;} // scale3 J/K
            real getElectricConversionFactor() override {return 8.987551e9*scale2;} // scale-2 N m^2 C^-2
    };

}}}

REGISTER_UNITS(
    Units,lenSI,
    uammd::structured::Units::lenSI
)
