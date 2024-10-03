Ensemble
========

If you want to implement a system where the macroscopic variables go further than
Volume and Temperature, you probably want to add a new ``Ensemble``. In practice, add
a new ``Ensemble`` is define new global variables that can be read and written from all
parts of the code. So you have to implement the private atributes that contains the global
variables, and the methods to write and read them. Also ``updateDataEntry()`` handles the
comunication with the dataEntry object.

.. code-block:: cpp

    #include "GlobalData/Ensemble/EnsembleHandler.cuh"
    #include "GlobalData/Ensemble/EnsembleFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Ensemble{

        class myEnsemble: public EnsembleHandler{

            private:

                Box box;
                real temperature;
                real myVariable;

            public:

                myEnsemble(DataEntry& data):EnsembleHandler(data){

                    auto ensembleData = data.getDataMap();

                    temperature   = ensembleData[0]["temperature"];

                    real3 boxSize = real3(ensembleData[0]["box"]);
                    box           = Box(boxSize);

                    myVariable        = ensembleData[0]["myVariable"];
                }

                Box  getBox()         override{return box;}
                real getTemperature() override{return temperature;}
                real getMyVariable()  override{return myVariable;}

                void setBox(Box newBox)                  override{box = newBox;}
                void setTemperature(real newTemperature) override{temperature = newTemperature;}
                void setmyVariable(real newMyVariable)   override{myVariable = newMyVariable;}

                void updateDataEntry(DataEntry data) override{
                    std::vector<std::string> labels = data.getLabels();
                    for(int i = 0; i < labels.size(); i++){
                        std::string lbl = labels[i];

                        if(lbl == "temperature"){
                            data.setData(0,i,this->temperature);
                        }

                        if(lbl == "box"){
                            data.setData(0,i,box.boxSize);
                        }

                        if(lbl == "myVariable"){
                            data.setData(0,i,this->myVariable);
                        }
                    }
                }

        };

    }}}

    REGISTER_ENSEMBLE(
        Ensemble,myEnsemble,
        uammd::structured::Ensemble::myEnsemble
    )

To register your own Ensemble system create the file
``src/GlobalData/Ensemble/myEnsemble/myEnsemble.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "GlobalData":
        "Ensemble":[
            ["..."],
            ["Ensemble","myEnsemble","myEnsemble.cu"]
            ]
   }

New Ensemble Global variables must be included in the ``Data.json`` file.

.. code-block:: json
   :emphasize-lines: 6

    "Ensemble": [
        ["lambda", "Lambda", "real"],
        ["temperature", "Temperature", "real"],
        ["box", "Box", "Box"],
        ["myVariable", "MyVariable", "DataType"]
    ]


