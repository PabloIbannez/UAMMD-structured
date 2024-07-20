#include "GlobalData/Ensemble/EnsembleHandler.cuh"
#include "GlobalData/Ensemble/EnsembleFactory.cuh"

namespace uammd{
namespace structured{
namespace Ensemble{

    class NVT: public EnsembleHandler{

        private:

            Box box;
            real temperature;

        public:

            NVT(DataEntry& data):EnsembleHandler(data){

                auto ensembleData = data.getDataMap();

                temperature   = ensembleData[0]["temperature"];

                real3 boxSize = real3(ensembleData[0]["box"]);
                box           = Box(boxSize);
            }

            Box  getBox()         override{return box;}
            real getTemperature() override{return temperature;}

            void setBox(Box newBox)                  override{box = newBox;}
            void setTemperature(real newTemperature) override{temperature = newTemperature;}

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
                }
            }

    };

}}}

REGISTER_ENSEMBLE(
    Ensemble,NVT,
    uammd::structured::Ensemble::NVT
)
