#include "GlobalData/Ensemble/EnsembleHandler.cuh"
#include "GlobalData/Ensemble/EnsembleFactory.cuh"

namespace uammd{
namespace structured{
namespace Ensemble{

    class NVTlambda: public EnsembleHandler{

        private:

            Box box;
            real temperature;
            real lambda;

        public:

            NVTlambda(DataEntry& data):EnsembleHandler(data){

                auto ensembleData = data.getDataMap();

                temperature   = ensembleData[0]["temperature"];

                real3 boxSize = real3(ensembleData[0]["box"]);
                box           = Box(boxSize);

                lambda        = ensembleData[0]["lambda"];
            }

            Box  getBox()         override{return box;}
            real getTemperature() override{return temperature;}
            real getLambda()      override{return lambda;}

            void setBox(Box newBox)                  override{box = newBox;}
            void setTemperature(real newTemperature) override{temperature = newTemperature;}
            void setLambda(real newLambda)           override{lambda = newLambda;}

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

                    if(lbl == "lambda"){
                        data.setData(0,i,this->lambda);
                    }
                }
            }

    };

}}}

REGISTER_ENSEMBLE(
    Ensemble,NVTlambda,
    uammd::structured::Ensemble::NVTlambda
)
