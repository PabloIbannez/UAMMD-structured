#ifndef __STRUCTURED__
#define __STRUCTURED__

#include <map>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <limits>

#include"ThirdParty/json.hpp"

#include"uammd.cuh"
#include"utils/container.h"

#include"misc/ParameterUpdatable.h"

#include"utils/InputFile.h"

#include"Interactor/ExternalForces.cuh"

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
//Virial -real(0.5)*dot(rij,Fij)
//Stress -real(0.5)*outer(rij,Fij)

namespace uammd{
namespace structured{
namespace quaternions{

    VECATTR real3 getEx(const real4& quat) { 
                                            real q0 = quat.x;
                                            real q1 = quat.y;
                                            real q2 = quat.z;
                                            real q3 = quat.w;

                                            return real(2.0)*make_real3(q0*q0+q1*q1-real(0.5),
                                                                        q1*q2+q0*q3,
                                                                        q1*q3-q0*q2);}

    VECATTR real3 getEy(const real4& quat) { 
                                            real q0 = quat.x;
                                            real q1 = quat.y;
                                            real q2 = quat.z;
                                            real q3 = quat.w;

                                            return real(2.0)*make_real3(q1*q2-q0*q3,
                                                                        q0*q0+q2*q2-real(0.5),
                                                                        q2*q3+q0*q1);}

    VECATTR real3 getEz(const real4& quat) { 
                                            real q0 = quat.x;
                                            real q1 = quat.y;
                                            real q2 = quat.z;
                                            real q3 = quat.w;

                                            return real(2.0)*make_real3(q1*q3+q0*q2,
                                                                        q2*q3-q0*q1,
                                                                        q0*q0+q3*q3-real(0.5));}

    VECATTR real4 rotate(const real4& quat,const real3& omega) {
                                                  
        real q0 = quat.x;
        real q1 = quat.y;
        real q2 = quat.z;
        real q3 = quat.w;

        real4 rot;

        real w1=omega.x;
        real w2=omega.y;
        real w3=omega.z;

        real factor = -dot(omega,omega)/real(8.0);

        rot.x = q0  +  real(0.5)*(-q1*w1-q2*w2-q3*w3)  +  factor*q0;
        rot.y = q1  +  real(0.5)*(+q0*w1+q3*w2-q2*w3)  +  factor*q1;
        rot.z = q2  +  real(0.5)*(-q3*w1+q0*w2+q1*w3)  +  factor*q2;
        rot.w = q3  +  real(0.5)*(+q2*w1-q1*w2+q0*w3)  +  factor*q3;

        return rot;

    }

}}}

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

__device__ __host__ real computeVirial(const real3& rij,const real3& Fij){
    return real(0.5)*dot(-rij,Fij);
}

__device__ __host__ tensor3 computeStress(const real3& rij,const real3& Fij){
    return real(0.5)*outer(-rij,Fij);
}

SFINAE_DEFINE_HAS_MEMBER(force);
SFINAE_DEFINE_HAS_MEMBER(energy);
SFINAE_DEFINE_HAS_MEMBER(virial);
SFINAE_DEFINE_HAS_MEMBER(stress);

SFINAE_DEFINE_HAS_MEMBER(getForceTransverser);
SFINAE_DEFINE_HAS_MEMBER(getEnergyTransverser);
SFINAE_DEFINE_HAS_MEMBER(getVirialTransverser);
SFINAE_DEFINE_HAS_MEMBER(getStressTransverser);

SFINAE_DEFINE_HAS_MEMBER(box);

SFINAE_DEFINE_HAS_MEMBER(isConstrained);
  
}}

#include"Utils/selectors.cuh"
#include"Utils/exclusionsList.cuh"
#include"Utils/conditions.cuh"
#include"Utils/groupUtils.cuh"

#include"Measures/MeasuresBasic.cuh"

#include"Potentials/Potentials.cuh"

#include"Interactor/BondedInteractor.cuh"
#include"Interactor/PairInteractor.cuh"
   
#include"Interactor/GenericSurfaceInteractor.cuh"
#include"Interactor/SphericalShellInteractor.cuh"

#include"Interactor/AFMInteractor.cuh"
#include"Interactor/PlatesInteractor.cuh"

#include"Interactor/SetsInteractor.cuh"

#include"Interactor/CenterOfMassInteractor.cuh"
#include"Interactor/UmbrellaInteractor.cuh"

#include"Constraint/Constraint.cuh"

#include"ForceFields/ForceFields.cuh"

#include"Integrator/IntegratorBasic.cuh"
#include"Integrator/steepestDescent.cuh"
#include"Integrator/langevinNVT.cuh"
#include"Integrator/brownianNVT.cuh"

#include"Simulation/SimulationStep.cuh"

#include"Simulation/Simulation.cuh"
#include"Simulation/SimulationAFM.cuh"
#include"Simulation/SimulationPlates.cuh"
#include"Simulation/SimulationSphere.cuh"
#include"Simulation/SimulationUmbrellaCenterOfMassDistance.cuh"
#include"Simulation/SimulationUmbrellaAlongVector.cuh"

#include"Wrapper/Wrapper.cuh"

#endif
