#ifndef __NVT_ENSEMBLE__
#define __NVT_ENSEMBLE__

namespace uammd{
namespace structured{
namespace Ensemble{

    class NVT: public EnsembleHandler{

        private:

            Box box;
            real temperature;

        public:

            NVT(DataEntry& data);

            Box  getBox()         override;
            real getTemperature() override;

            void setBox(Box newBox)                  override;
            void setTemperature(real newTemperature) override;

            void updateDataEntry(DataEntry data) override;

    };

}}}

#endif
