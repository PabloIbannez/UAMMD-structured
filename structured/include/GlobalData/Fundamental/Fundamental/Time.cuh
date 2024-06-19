#ifndef __TIME_FUNDAMENTAL__
#define __TIME_FUNDAMENTAL__

namespace uammd{
namespace structured{
namespace Fundamental{

    class Time : public FundamentalHandler {

    private:

        ullint currentStep;
        double timeStep;
        double simulationTime;

        int nTimeSteps = 0;

    public:
        Time(DataEntry& data);

        ullint getCurrentStep() override;
        double getTimeStep() override;
        double getSimulationTime() override;

        void setCurrentStep(ullint newStep) override;
        void setTimeStep(double newTimeStep) override;
        void setSimulationTime(double newSimulationTime) override;

        void updateDataEntry(DataEntry data);
    };

}}}

#endif
