#ifndef __DYNAMICALLY_BONDED_PATCHY_PARTICLES_FUNDAMENTAL__
#define __DYNAMICALLY_BONDED_PATCHY_PARTICLES_FUNDAMENTAL__

namespace uammd {
namespace structured {
namespace Fundamental {

class DynamicallyBondedPatchyParticles : public FundamentalHandler {

private:

    double energyThreshold;

public:
    DynamicallyBondedPatchyParticles(DataEntry& data);

    void setEnergyThreshold(real newEnergyThreshold) override;
    real getEnergyThreshold() override;

    void updateDataEntry(DataEntry data);
};

}}}

#endif
