#include "GlobalData/Ensemble/EnsembleHandler.cuh"
#include "GlobalData/Ensemble/Ensemble/NVTlambda.cuh"

namespace uammd{
namespace structured{
namespace Ensemble{

    NVTlambda::NVTlambda(DataEntry& data):EnsembleHandler(data){

        auto ensembleData = data.getDataMap();

        temperature   = ensembleData[0]["temperature"];

        real3 boxSize = real3(ensembleData[0]["box"]);
        box           = Box(boxSize);

        lambda        = ensembleData[0]["lambda"];
    }

    Box NVTlambda::getBox()          { return box; }
    real NVTlambda::getTemperature() { return temperature; }
    real NVTlambda::getLambda()      { return lambda; }

    void NVTlambda::setBox(Box newBox)                  { box = newBox; }
    void NVTlambda::setTemperature(real newTemperature) { temperature = newTemperature; }
    void NVTlambda::setLambda(real newLambda)           { lambda = newLambda; }

    void NVTlambda::updateDataEntry(DataEntry data) {
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

}}}
