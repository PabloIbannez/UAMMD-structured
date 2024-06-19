#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/Fundamental/Time.cuh"

namespace uammd {
namespace structured {
namespace Fundamental {

    Time::Time(DataEntry& data) : FundamentalHandler(data) {
        this->setCurrentStep(data.getParameter<ullint>("currentStep", 0));
        this->setSimulationTime(data.getParameter<double>("simulationTime", 0.0));
    
        if (data.isParameterAdded("timeStep")) {
            this->setTimeStep(data.getParameter<double>("timeStep"));
        } else {
            System::log<System::WARNING>("[Time] No timeStep specified, using 0.0 as default.");
            timeStep = 0.0;
        }
    }
    
    ullint Time::getCurrentStep() {
        return currentStep;
    }
    
    double Time::getTimeStep() {
        if (nTimeSteps > 1) {
            System::log<System::CRITICAL>("[Time] Time step has been changed more than once. Not should be used");
        }
        return timeStep;
    }
    
    double Time::getSimulationTime() {
        return simulationTime;
    }
    
    void Time::setCurrentStep(ullint newStep) {
        currentStep = newStep;
    }
    
    void Time::setTimeStep(double newTimeStep) {
        timeStep = newTimeStep;
        nTimeSteps += 1;
    }
    
    void Time::setSimulationTime(double newSimulationTime) {
        simulationTime = newSimulationTime;
    }
    
    void Time::updateDataEntry(DataEntry data) {
        data.setParameter("currentStep", currentStep);
        data.setParameter("timeStep", timeStep);
        data.setParameter("simulationTime", simulationTime);
    }

}}}






