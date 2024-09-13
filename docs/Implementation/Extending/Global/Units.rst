Units
=====

To set your own units system you must define physical constants
in your units. Boltzman constant :math:`k_bT` and electric conversion
factor :math:`\frac{1}{4\pi\epsilon_0}` are the current constant you must define.

.. code-block:: cpp

    #include "GlobalData/Units/UnitsHandler.cuh"
    #include "GlobalData/Units/UnitsFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Units{

        class myUnits: public UnitsHandler{

            public:

                None(DataEntry& data):UnitsHandler(data){}

                real getBoltzmannConstant()        override {return myUnits_KbT;}
                real getElectricConversionFactor() override {return myUnits_Ke;}
        };

    }}}

    REGISTER_UNITS(
        Units,myUnits,
        uammd::structured::Units::myUnits
    )

To register your own Units system create the file
``src/GlobalData/Units/Units/myUnits.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "GlobalData":
        "Units":[
            ["..."],
            ["Units","myUnits","myUnits.cu"]
            ]
   }

If you want to add a new constant, you must first add it to the "Units" section of ``structured/Data.json``,
and then just implement ``getMyConstant()`` in the ``myUnits`` class. If
you do that, take into account that every other Units system will raise
an error when you ask for MyConstant, since they don't have ``getMyConstant()`` defined.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Units": [
       ["BoltzmannConstant", "BoltzmannConstant", "real"],
       ["electricConversionFactor", "ElectricConversionFactor", "real"],
       ["MyConstant", "MyConstant","DataType"]
   }

