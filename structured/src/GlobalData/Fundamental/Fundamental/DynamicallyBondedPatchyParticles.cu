#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/Fundamental/DynamicallyBondedPatchyParticles.cuh"

namespace uammd {
namespace structured {
namespace Fundamental {

    DynamicallyBondedPatchyParticles::DynamicallyBondedPatchyParticles(DataEntry& data) : FundamentalHandler(data) {
        this->setEnergyThreshold(data.getParameter<real>("energyThreshold", 0.0));
    }
    
    void DynamicallyBondedPatchyParticles::setEnergyThreshold(real newEnergyThreshold) {
        energyThreshold = newEnergyThreshold;
    }
    
    real DynamicallyBondedPatchyParticles::getEnergyThreshold() {
        return energyThreshold;
    }
    
    void DynamicallyBondedPatchyParticles::updateDataEntry(DataEntry data) {
        data.setParameter("energyThreshold", energyThreshold);
    }

}}}
