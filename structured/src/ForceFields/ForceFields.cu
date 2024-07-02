#include "ForceFields/ForceFields.cuh"

namespace uammd{
namespace structured{
namespace ForceField_ns{

    void Scheduled::stopInteractor(std::string interactorName){

        ullint currentStep = gd->getFundamental()->getCurrentStep();

        auto interactor = interactors.find(interactorName);
        if(interactor == interactors.end()){
            System::log<System::CRITICAL>("[ForceField] (%s) Trying to stop an interactor that does not exist: %s",
                                          this->name.c_str(),interactorName.c_str());
        }
        idleInteractors[interactorName] = interactor->second;
        interactors.erase(interactor);

        System::log<System::MESSAGE>("[ForceField] (%s) Stopped interactor: %s, at step: %llu",
                                      this->name.c_str(),interactorName.c_str(),currentStep);
    }

    void Scheduled::resumeInteractor(std::string interactorName){

        ullint currentStep = gd->getFundamental()->getCurrentStep();

        auto interactor = idleInteractors.find(interactorName);
        if(interactor == idleInteractors.end()){
            System::log<System::CRITICAL>("[ForceField] (%s) Trying to resume an interactor that has not been stopped: %s",
                                          this->name.c_str(),interactorName.c_str());
        }
        interactors[interactorName] = interactor->second;
        idleInteractors.erase(interactor);

        System::log<System::MESSAGE>("[ForceField] (%s) Resumed interactor: %s, at step: %llu",
                                     this->name.c_str(),interactorName.c_str(),currentStep);
    }

    Scheduled::Scheduled(std::shared_ptr<Topology>   top):Base(top,"Scheduled"){

        gd = top->getGlobalData();

        for(auto &interactor: interactors){
            DataEntry data = top->getSystem()->getInput()->getDataEntry(top->getInteractorInfo(interactor.first).path);
            if(data.isParameterAdded("endStep") or data.isParameterAdded("startStep")){
                schedule sched;
                sched.start = data.getParameter<ullint>("startStep",0);
                sched.end   = data.getParameter<ullint>("endStep",std::numeric_limits<ullint>::max());

                if(sched.start > sched.end){
                    System::log<System::CRITICAL>("[ForceField] (%s) Interactor %s has startStep > endStep",
                                                  this->name.c_str(),interactor.first.c_str());
                }

                scheduledInteractors[interactor.first] = sched;
            }
        }

        ullint step = this->gd->getFundamental()->getCurrentStep();
        for(auto &interactor: scheduledInteractors){
            if(interactor.second.start <= step and step <= interactor.second.end){
                interactor.second.state = true;
            }
            else{
                interactor.second.state = false;
                this->stopInteractor(interactor.first);
            }
        }
    }

    void Scheduled::sum(Computables comp,cudaStream_t st) {
        ullint step = this->gd->getFundamental()->getCurrentStep();
        for(auto &interactor: scheduledInteractors){
            if(interactor.second.start <= step and step < interactor.second.end){
                if(not interactor.second.state){
                    this->resumeInteractor(interactor.first);
                    interactor.second.state = true;
                }
            }
            else{
                if(interactor.second.state){
                    this->stopInteractor(interactor.first);
                    interactor.second.state = false;
                }
            }
        }

        Base::sum(comp,st);
    }

}}}
