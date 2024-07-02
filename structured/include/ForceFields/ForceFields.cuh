#pragma once

#include "Topology/Topology.cuh"
#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

namespace uammd{
namespace structured{

namespace ForceField_ns{

using Computables = uammd::Interactor::Computables;

class ForceFieldBase: public Interactor{

        std::shared_ptr<Topology>    top;

    public:

        ForceFieldBase(std::shared_ptr<Topology>   top,
                       std::string name):Interactor(top->getParticleGroup(),
                                                    name),top(top){}

        std::shared_ptr<Topology>   getTopology(){return top;}

};

class None : public ForceFieldBase{

    private:

        using Base = ForceFieldBase;

    public:

        None(std::shared_ptr<Topology>   top):Base(top,"None"){}

        void sum(Computables comp,cudaStream_t st) override {return;}
};

class Generic : public ForceFieldBase {

    protected:

        using Base = ForceFieldBase;

        std::map<std::string, std::shared_ptr<Interactor>> interactors;
        std::map<std::string, std::shared_ptr<Interactor>> idleInteractors;

    public:

        Generic(std::shared_ptr<Topology>   top,std::string name):Base(top,name){
            interactors = top->getInteractors();
        }

        Generic(std::shared_ptr<Topology>   top):Generic(top,"Generic"){}

        void sum(Computables comp,cudaStream_t st) override {
            for(auto &interactor: interactors){
                interactor.second->sum(comp,st);
            }
        }

};

class Scheduled : public Generic {

    private:

        struct schedule{
            ullint start;
            ullint end;
            bool state; //true: on, false: off
        };

        using Base = Generic;

        std::shared_ptr<GlobalData> gd;

        std::map<std::string,schedule> scheduledInteractors;

        void stopInteractor(std::string interactorName);
        void resumeInteractor(std::string interactorName);

    public:

        Scheduled(std::shared_ptr<Topology>   top);

        void sum(Computables comp,cudaStream_t st);

};

}

using ForceField = ForceField_ns::Scheduled;

}}
