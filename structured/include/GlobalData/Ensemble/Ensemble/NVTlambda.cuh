#ifndef __NVTLAMBDA_ENSEMBLE__
#define __NVTLAMBDA_ENSEMBLE__

namespace uammd{
namespace structured{
namespace Ensemble{

    class NVTlambda: public EnsembleHandler{

        private:

            Box box;
            real temperature;
            real lambda;

        public:

            NVTlambda(DataEntry& data);

            Box  getBox()         override;
            real getTemperature() override;
            real getLambda()      override;

            void setBox(Box newBox)                  override;
            void setTemperature(real newTemperature) override;
            void setLambda(real newLambda)           override;

            void updateDataEntry(DataEntry data) override;

    };

}}}

#endif

