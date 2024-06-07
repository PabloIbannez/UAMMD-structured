#ifndef __DYNAMICALLY_BONDED_PATCHY_PARTICLES_FUNDAMENTAL__
#define __DYNAMICALLY_BONDED_PATCHY_PARTICLES_FUNDAMENTAL__

namespace uammd{
namespace structured{
namespace Fundamental{

    class DynamicallyBondedPatchyParticles: public FundamentalHandler{

        private:

            double energyThreshold;

        public:

            DynamicallyBondedPatchyParticles(DataEntry& data):FundamentalHandler(data){

                this->setEnergyThreshold(data.getParameter<real>("energyThreshold",0.0));
            }

            void setEnergyThreshold(real newEnergyThreshold) override{energyThreshold = newEnergyThreshold;}

            real getEnergyThreshold()  override {return energyThreshold;}

            void updateDataEntry(DataEntry data){
                data.setParameter("energyThreshold",energyThreshold);
            }

    };

}}}

#endif
