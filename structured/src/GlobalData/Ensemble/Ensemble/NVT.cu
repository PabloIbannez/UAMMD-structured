#include "GlobalData/Ensemble/EnsembleHandler.cuh"
#include "GlobalData/Ensemble/Ensemble/NVT.cuh"

namespace uammd{
namespace structured{
namespace Ensemble{

    NVT::NVT(DataEntry& data):EnsembleHandler(data){

        auto ensembleData = data.getDataMap();

        temperature   = ensembleData[0]["temperature"];

        real3 boxSize = real3(ensembleData[0]["box"]);
        box           = Box(boxSize);
    }

    Box NVT::getBox()          { return box; }
    real NVT::getTemperature() { return temperature; }

    void NVT::setBox(Box newBox)                  { box = newBox; }
    void NVT::setTemperature(real newTemperature) { temperature = newTemperature; }

    void NVT::updateDataEntry(DataEntry data) {
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

}}}

