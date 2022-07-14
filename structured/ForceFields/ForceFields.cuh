#ifndef __FORCE_FIELD__
#define __FORCE_FIELD__

namespace uammd{
namespace structured{
namespace forceField{

using Computables = uammd::Interactor::Computables;

template<class Units_,
         class Types_ >
class ForceFieldBase: public Interactor{

    public:

        using Units = Units_;
        using Types = Types_;

        using Topology  = Topology<Units_,Types_>;

    protected:

        std::shared_ptr<Topology> top;

    public:

        ForceFieldBase(std::shared_ptr<ParticleGroup> pg,
                       InputFile&                     in):Interactor(pg,"ForceField"){
            top = std::make_shared<Topology>(sys,in);                                     
        }
        
        std::shared_ptr<Topology> getTopology(){return top;};

        void sum(Computables comp,cudaStream_t st) override {return;}


};

template<class Units_,
         class Types_,
         template <class Topology_> class Condition_>
class ForceFieldNeighbourBase : public ForceFieldBase<Units_,Types_>{

    protected:

        using Base = ForceFieldBase<Units_,Types_>;

        using Condition = Condition_<typename Base::Topology>;
        
        using NeighbourList = ConditionedVerletListSet<Condition>;
        
        std::shared_ptr<Condition> condition;
        std::shared_ptr<NeighbourList>    nl; 
            
        real VerletListDst;

    public:

        ForceFieldNeighbourBase(std::shared_ptr<ParticleGroup> pg,
                                InputFile&                     in):Base(pg,in),
                                                                   VerletListDst(std::stof(in.getOption("VerletListDst",InputFile::Required).str())){
            
            this->sys->template log<System::MESSAGE>("[ForceFieldNeighbourBase] "
                                                     "Parameter VerletListDst added: %f",
                                                      VerletListDst);
                

            condition = std::make_shared<Condition>(this->pd,this->top,in);
            
            typename NeighbourList::Parameters NeighbourListParam;
            
            NeighbourListParam.cutOff       = real(0.0);
            NeighbourListParam.cutOffVerlet = VerletListDst;

            nl = std::make_shared<NeighbourList>(this->pg,
                                                 condition,
                                                 NeighbourListParam);
        
        
        }

        void sum(Computables comp,cudaStream_t st) override {return;}

};


template<class Units_,
         class Types_ >
class none : public ForceFieldBase<Units_,Types_>{

    private:

        using Base = ForceFieldBase<Units_,Types_>;

        std::vector<std::string> componentsList;

    public:

        none(std::shared_ptr<ParticleGroup> pg,
             InputFile&                     in):Base(pg,in){}
        
        std::vector<std::string> getComponentsList(){return componentsList;}

        void sum(Computables comp,cudaStream_t st) override {return;}
        
        void sum(std::string component,Computables comp,cudaStream_t st) {
            this->sys->template log<System::CRITICAL>("[NONE] Requested potential %s to sum. "
                                                        "But %s is not present in the force field",
                                                        component.c_str(),component.c_str());
        }
};

}}}

#include"Generic/Generic.cuh"

#endif
