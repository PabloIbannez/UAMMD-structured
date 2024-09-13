Fundamental
===========

If you have already read the section on ``Fundamental``, it is hard to imagine, 
even for the writer, what has led you to open this section. Adding a new 
Fundamental means that you are going to do something very different with 
UAMMD-structured than it was intended for. You can, for example, convert
UAMMD-structured in a Variational Method solver, so there your Fundamental
variable is not **Time** but the **Iterative Error**.

.. warning::
    Many of the UAMMD-structured 
    modules may NOT work if time is not the Fundamental variable in your simulation. 
    For example, most of the Integrators (if not all) request the time from 
    GlobalData to add dt and thus integrate.

Only mandatory function of ``Fundamental`` is ``updateDataEntry(DataEntry data)``
that contains the logic of how to update the fundamental variable(s). But Fundamental
class should be decorated with private atributes (fundamental variables) and the public
methods to handle them. In the example only one fundamental variable ``myVariable`` is
implemented. You will have access to the class Fundamental in other parts of the code by
doing ``gd->getFundamental()``. So your can read and write your Fundamental variable(s).

.. code-block:: cpp

    #include "GlobalData/Fundamental/FundamentalHandler.cuh"
    #include "GlobalData/Fundamental/FundamentalFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Fundamental{

        class myFundamental: public FundamentalHandler{
            private:
                double myVariable;

            public:

                myFundamental(DataEntry& data):FundamentalHandler(data){
                    this-setMyVariable(data.getParameter<real>("myVariable",0.0)); //if myVariable not in input myVariable is set to 0.0
                }

                setMyVariable(real newMyVariable) override{myVariable = newMyVariable;}

                getMyVariable()                   override{return myVariable;}

                void updateDataEntry(DataEntry data){
                    data.setParameter("myVariable", myVariable);
                }

        };

    }}}

    REGISTER_FUNDAMENTAL(
        Fundamental,myFundamental,
        uammd::structured::Fundamental::myFundamental
    )
