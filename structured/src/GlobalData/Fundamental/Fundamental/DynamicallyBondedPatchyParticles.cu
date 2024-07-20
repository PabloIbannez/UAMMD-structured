#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/FundamentalFactory.cuh"

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

REGISTER_FUNDAMENTAL(
    Fundamental,DynamicallyBondedPatchyParticles,
    uammd::structured::Fundamental::DynamicallyBondedPatchyParticles
)
