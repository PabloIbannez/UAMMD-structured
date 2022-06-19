#ifndef __WRAPPER__
#define __WRAPPER__

#include <map>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <limits>

#include"ThirdParty/json.hpp"

#include"uammd.cuh"

#include"misc/ParameterUpdatable.h"

#include"utils/InputFile.h"

#include"Interactor/ExternalForces.cuh"

#include"Utils/tensor.cuh"

#include"Utils/misc.cuh"
#include"Utils/parameterHandler.cuh"

#include"InputOutput/psf.cuh"
#include"InputOutput/dcd.cuh"

#include"InputOutput/InputOutput.cuh"

#include"Topology/Topology.cuh"
#include"Types/Types.cuh"

#include"NeighbourList/ConditionedListSetBase.cuh"
#include"NeighbourList/ConditionedVerletListSet.cuh"

#include"NeighbourList/FixedNeighbourList.cuh"

//GROMACS convention
//rij = rj - ri
//Fij force over particle i due to j
//Virial -real(0.5)*outer(rij,Fij)

namespace uammd{
namespace structured{
namespace UnitsSystem{

    struct NONE{
        
        inline static const std::string NAME = "NONE";
        
        inline static const std::string L_INTERNAL = "NONE";
        inline static const std::string L_EXTERNAL = "NONE";
        inline static const std::string T_INTERNAL = "NONE";
        inline static const std::string T_EXTERNAL = "NONE";
        inline static const std::string E_INTERNAL = "NONE";
        inline static const std::string E_EXTERNAL = "NONE";
        inline static const std::string F_INTERNAL = "NONE";
        inline static const std::string F_EXTERNAL = "NONE";
        
        static constexpr real FROM_INTERNAL_TIME = real(1.0);
        static constexpr real TO_INTERNAL_TIME   = real(1.0);
        
        static constexpr real TO_INTERNAL_FORCE   = real(1.0);
        static constexpr real FROM_INTERNAL_FORCE = real(1.0);
        
        static constexpr real KBOLTZ  = real(1.0);
        static constexpr real ELECOEF = real(1.0);

    };

    struct KCALMOL_A{
        //External kcal A piconewton picosecond
        //Internal AKMA

        inline static const std::string NAME = "KCALMOL_A";
        
        inline static const std::string L_INTERNAL = "A";
        inline static const std::string L_EXTERNAL = "A";
        inline static const std::string T_INTERNAL = "AKMA";
        inline static const std::string T_EXTERNAL = "ps";
        inline static const std::string E_INTERNAL = "Kcal/mol";
        inline static const std::string E_EXTERNAL = "Kcal/mol";
        inline static const std::string F_INTERNAL = "(Kcal/mol)/A";
        inline static const std::string F_EXTERNAL = "pN";
    
        static constexpr real TO_INTERNAL_TIME   = real(1.0)/real(4.88882129E-02);
        static constexpr real FROM_INTERNAL_TIME = real(1.0)/TO_INTERNAL_TIME;
        
        static constexpr real TO_INTERNAL_FORCE   = real(1.0)/real(69.5118);
        static constexpr real FROM_INTERNAL_FORCE = real(1.0)/TO_INTERNAL_FORCE;
        
        static constexpr real KBOLTZ  = 1.987191E-03;
        static constexpr real ELECOEF = 332.0716;
    };
}

__device__ __host__ tensor3 computeVirial(const real3& rij,const real3& Fij){
    return real(0.5)*outer(-rij,Fij);
}

SFINAE_DEFINE_HAS_MEMBER(box);
SFINAE_DEFINE_HAS_MEMBER(isConstrained);
  
}}

#include"Utils/selectors.cuh"
#include"Utils/exclusionsList.cuh"
#include"Utils/conditions.cuh"
#include"Utils/groupUtils.cuh"

//#include"Measures/MeasuresBasic.cuh"

#include"Potentials/Potentials.cuh"

#include"Interactor/BondedInteractor.cuh"
#include"Interactor/PairInteractor.cuh"

//#include"Interactor/AFMInteractor.cuh"
//#include"Interactor/QCMInteractor.cuh"
//#include"Interactor/PlatesInteractor.cuh"
//
//#include"Interactor/SetsInteractor.cuh"
//
//#include"Interactor/CenterOfMassInteractor.cuh"
//#include"Interactor/UmbrellaInteractor.cuh"
//
//#include"Interactor/GAMD.cuh"

#include"Constraint/Constraint.cuh"

#include"ForceFields/ForceFields.cuh"

namespace uammd{
namespace structured{
namespace Wrapper{

std::shared_ptr<ParticleData> setUpParticleData(std::shared_ptr<System> sys,
                                                InputFile& in){
    
    std::string inputCoordPath;

    in.getOption("inputCoordPath"   ,InputFile::Required) >> inputCoordPath;

    struct coordFormat{
        int   id;
        real3 pos;
        real3 vel;
        real4 dir;
    };
    
    std::vector<coordFormat> pdBuffer = InputOutput::loadCoordFromFile<coordFormat>(sys,inputCoordPath);
            
    std::shared_ptr<ParticleData> pd = std::make_shared<ParticleData>(pdBuffer.size(),sys);
            
    auto pId = pd->getId(access::location::cpu, access::mode::write);
    auto pos = pd->getPos(access::location::cpu, access::mode::write);
    auto vel = pd->getVel(access::location::cpu, access::mode::write);
    auto dir = pd->getDir(access::location::cpu, access::mode::write);

    fori(0,pdBuffer.size()){

        if(i==pdBuffer[i].id){
            pId[i]   = pdBuffer[i].id;
            pos[i].x = pdBuffer[i].pos.x;
            pos[i].y = pdBuffer[i].pos.y;
            pos[i].z = pdBuffer[i].pos.z;
            vel[i].x = pdBuffer[i].vel.x;
            vel[i].y = pdBuffer[i].vel.y;
            vel[i].z = pdBuffer[i].vel.z;
            dir[i].x = pdBuffer[i].dir.x;
            dir[i].y = pdBuffer[i].dir.y;
            dir[i].z = pdBuffer[i].dir.z;
            dir[i].w = pdBuffer[i].dir.w;

        } else {
            sys->log<System::CRITICAL>("[Wrapper] The internal id has to "
                    "match with the given by coord file. Inte: %i, File %i",i,pdBuffer[i].id);
        }
    }

    return pd;
}

std::shared_ptr<ParticleGroup> setUpParticleGroup(std::shared_ptr<System>       sys,
                                                  std::shared_ptr<ParticleData>  pd,
                                                  InputFile& in){
            
    std::shared_ptr<ParticleGroup> pg = std::make_shared<ParticleGroup>(pd,"All");

    return pg;
}


template<class ForceField>        
std::shared_ptr<ForceField> setUpForceField(std::shared_ptr<System>       sys,
                                            std::shared_ptr<ParticleData>  pd,
                                            std::shared_ptr<ParticleGroup> pg,
                                            InputFile& in){
    
    std::shared_ptr<ForceField> ff = std::make_shared<ForceField>(sys,pd,pg,in);
            
    auto top = ff->getTopology();

    top->loadStructureData(pd);
    top->loadTypes(pd);

    return ff;
}

}}}

#endif
