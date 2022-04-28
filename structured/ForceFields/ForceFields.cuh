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

        ForceFieldBase(std::shared_ptr<System>        sys,
                       std::shared_ptr<ParticleData>  pd,
                       std::shared_ptr<ParticleGroup> pg,
                       InputFile&                     in):Interactor(pd,pg,sys,"ForceField"){
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

        ForceFieldNeighbourBase(std::shared_ptr<System>        sys,
                                std::shared_ptr<ParticleData>  pd,
                                std::shared_ptr<ParticleGroup> pg,
                                InputFile&                     in):Base(sys,pd,pg,in),
                                                                   VerletListDst(std::stof(in.getOption("VerletListDst",InputFile::Required).str())){
            
            this->sys->template log<System::MESSAGE>("[ForceFieldNeighbourBase] "
                                                     "Parameter VerletListDst added: %f",
                                                      VerletListDst);
                

            condition = std::make_shared<Condition>(this->sys,this->pd,this->top,in);
            
            typename NeighbourList::Parameters NeighbourListParam;
            
            NeighbourListParam.cutOff       = real(0.0);
            NeighbourListParam.cutOffVerlet = VerletListDst;

            nl = std::make_shared<NeighbourList>(this->sys,this->pd,this->pg,
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

        none(std::shared_ptr<System>        sys,
             std::shared_ptr<ParticleData>  pd,
             std::shared_ptr<ParticleGroup> pg,
             InputFile&                     in):Base(sys,pd,pg,in){}
        
        std::vector<std::string> getComponentsList(){return componentsList;}

        void sum(Computables comp,cudaStream_t st) override {return;}
        
        void sum(std::string component,Computables comp,cudaStream_t st) {
            this->sys->template log<System::CRITICAL>("[NONE] Requested potential %s to sum. "
                                                        "But %s is not present in the force field",
                                                        component.c_str(),component.c_str());
        }
};

}}}

#include"ElasticNetworkModel/ElasticNetworkModel.cuh"

#include"Disordered/LinkerBonded.cuh"

#include"WeeksChandlerAndersen/WeeksChandlerAndersen.cuh"
#include"WeeksChandlerAndersen/WeeksChandlerAndersenLinkers.cuh"

#include"LennardJones/LennardJones.cuh"
#include"LJ_WCA/LJ_WCA.cuh"

#include"Generic/GenericBonded.cuh"

#include"Statistical/Statistical.cuh"

#include"Polymer/KratkyPorod.cuh"
#include"Polymer/Polymer.cuh"

#include"ShapeBasedCoarseGrained/ShapeBasedCoarseGrained.cuh"

#include"Clashed/Clashed.cuh"
#include"Clashed/ClashedLinkers.cuh"

#include"KaranicolasBrooks/KaranicolasBrooks.cuh"
#include"SelfOrganizedPolymer/SelfOrganizedPolymer.cuh"
#include"MechanicallyAccurateDNA/MechanicallyAccurateDNA.cuh"

#include"KimHummer/KimHummer.cuh"
#include"KimHummer/KimHummerLinkers.cuh"

#include"Ravikumar/Ravikumar.cuh"

#include"Bounds/SphericalShell.cuh"
#include"LennardJonesSurface/LennardJonesSurface.cuh"
#include"WeeksChandlerAndersenSurface/WeeksChandlerAndersenSurface.cuh"

#include"KaranicolasBrooks/KaranicolasBrooksSurface.cuh"

//#include"DLVO/DLVO.cuh"

namespace uammd{
namespace structured{
namespace forceField{

    //
    using EMPTY_KCALMOL_A = none<UnitsSystem::KCALMOL_A,
                                 Types::BASIC>;
    
    //WCA
    using WCA = WeeksChandlerAndersen::WeeksChandlerAndersen<ForceFieldNeighbourBase<UnitsSystem::NONE,
                                                                                     Types::BASIC,
                                                                                     conditions::allInter>>; 

    //LJ
    using LJ  = LennardJones::LennardJones<ForceFieldNeighbourBase<UnitsSystem::NONE,
                                                                   Types::BASIC,
                                                                   conditions::allInter>>; 
    
    //SBCG
    using SBCG = ShapeBasedCoarseGrained::ShapeBasedCoarseGrained<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                          Types::BASIC,
                                                                                          conditions::excluded>>;

    //KP
    using KP   = KratkyPorodModel::KratkyPorodModel<ForceFieldBase<UnitsSystem::NONE,Types::BASIC>>;

    //GenericBonded
    using GB_BASE = Generic::GenericBonded<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                   Types::BASIC,
                                                                   conditions::excludedIntraInterCharged>>;

    //Statistical
    using STAT    = Statistical::Statistical<GB_BASE>;

    //MechanicallyAccurateDNA
    using MADna = MechanicallyAccurateDNA::MechanicallyAccurateDNA<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                           Types::BASIC,
                                                                                           conditions::chargedExcluded>>;
    
    //KaranicolasBrooks
    using KB = KaranicolasBrooks::KaranicolasBrooks<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                            Types::BASIC,
                                                                            conditions::excludedIntra>>;
    
    //KaranicolasBrooksSurface
    using KB_SURF = KaranicolasBrooks::KaranicolasBrooksSurface<KaranicolasBrooks::KaranicolasBrooks<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                                                             Types::SURF,
                                                                                                                             conditions::excludedIntra>>>;
    
    //SelfOrganizedPolymer
    using SOP       = SelfOrganizedPolymer::SelfOrganizedPolymer<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                         Types::BASIC,
                                                                                         conditions::excludedIntra>>;
    
    using SOPG     = SelfOrganizedPolymer::SelfOrganizedPolymerGaussian<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                                Types::BASIC,
                                                                                                conditions::excludedIntra>>;

    using SOP_FIXED = SelfOrganizedPolymer::SelfOrganizedPolymerFixedSteric<ForceFieldBase<UnitsSystem::KCALMOL_A,
                                                                                           Types::BASIC>>;
    
    using SOP_BASE = SelfOrganizedPolymer::SelfOrganizedPolymer<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                        Types::BASIC,
                                                                                        conditions::excludedIntraInter>>;

    //Elastic Network Model
    using ENM_BASE = ElasticNetworkModel::ElasticNetworkModel<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                      Types::BASIC,
                                                                                      conditions::excludedIntraInter>>;
    
    using ENM_WCA   = WeeksChandlerAndersen::WeeksChandlerAndersen<ENM_BASE>;
    
    //Clashed
    using SOP_CLS   = Clashed::Clashed<SOP_BASE>;

    using SOP_CLS_L   = Clashed::ClashedLinkers<SOP_BASE>;

    using ENM_CLS   = Clashed::Clashed<ENM_BASE>;
    
    //

    using ENM_SOP_WCA   = WeeksChandlerAndersen::WeeksChandlerAndersen<SelfOrganizedPolymer::SelfOrganizedPolymer<ENM_BASE>>;

    //KimHummer
    
    using KHtype_A = Potentials::UnBound::KimHummerPotential_ns::ModelType::A;
    using KHtype_B = Potentials::UnBound::KimHummerPotential_ns::ModelType::B;
    using KHtype_C = Potentials::UnBound::KimHummerPotential_ns::ModelType::C;
    using KHtype_D = Potentials::UnBound::KimHummerPotential_ns::ModelType::D;
    using KHtype_E = Potentials::UnBound::KimHummerPotential_ns::ModelType::E;
    using KHtype_F = Potentials::UnBound::KimHummerPotential_ns::ModelType::F;

    using KHtype_LQLHK = Potentials::UnBound::KimHummerPotential_ns::ModelType::LQLHK;
    
    using SOP_KH_BASE = SelfOrganizedPolymer::SelfOrganizedPolymer<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                           Types::SASA,
                                                                                           conditions::excludedIntraInter>>;
    
    using ENM_KH_BASE = ElasticNetworkModel::ElasticNetworkModel<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                         Types::SASA,
                                                                                         conditions::excludedIntraInter>>;

    using SOP_KH_A = KimHummer::KimHummer<SOP_KH_BASE,KHtype_A>;
    using SOP_KH_B = KimHummer::KimHummer<SOP_KH_BASE,KHtype_B>;
    using SOP_KH_C = KimHummer::KimHummer<SOP_KH_BASE,KHtype_C>;
    using SOP_KH_D = KimHummer::KimHummer<SOP_KH_BASE,KHtype_D>;
    using SOP_KH_E = KimHummer::KimHummer<SOP_KH_BASE,KHtype_E>;
    using SOP_KH_F = KimHummer::KimHummer<SOP_KH_BASE,KHtype_F>;
    
    using SOP_KH_LQLHK = KimHummer::KimHummer<SOP_KH_BASE,KHtype_LQLHK>;
    
    using ENM_SOP_KH_A = KimHummer::KimHummer<SelfOrganizedPolymer::SelfOrganizedPolymer<ENM_KH_BASE>,KHtype_A>;
    
    using SOP_KH_A_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_A>;
    using SOP_KH_B_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_B>;
    using SOP_KH_C_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_C>;
    using SOP_KH_D_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_D>;
    using SOP_KH_E_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_E>;
    using SOP_KH_F_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_F>;
    
    using SOP_KH_LQLHK_L = KimHummer::KimHummerLinkers<SOP_KH_BASE,KHtype_LQLHK>;
    
    using ENM_KH_BASE = ElasticNetworkModel::ElasticNetworkModel<ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                                                                         Types::SASA,
                                                                                         conditions::excludedIntraInter>>;
    
    using ENM_KH_A = KimHummer::KimHummer<ENM_KH_BASE,KHtype_A>;
    using ENM_KH_B = KimHummer::KimHummer<ENM_KH_BASE,KHtype_B>;
    using ENM_KH_C = KimHummer::KimHummer<ENM_KH_BASE,KHtype_C>;
    using ENM_KH_D = KimHummer::KimHummer<ENM_KH_BASE,KHtype_D>;
    using ENM_KH_E = KimHummer::KimHummer<ENM_KH_BASE,KHtype_E>;
    using ENM_KH_F = KimHummer::KimHummer<ENM_KH_BASE,KHtype_F>;
    
    using ENM_KH_LQLHK = KimHummer::KimHummer<ENM_KH_BASE,KHtype_LQLHK>;
    
    //SOPWCA

    using SOP_WCA   = WeeksChandlerAndersen::WeeksChandlerAndersen<SOP_BASE>;
    using SOP_WCA_L = WeeksChandlerAndersen::WeeksChandlerAndersenLinkers<SOP_BASE>;
    

}}}


#endif
