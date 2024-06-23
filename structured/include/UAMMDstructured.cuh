#ifndef __STRUCTURED__
#define __STRUCTURED__

#include"UAMMDstructuredBase.cuh"

//Define UAMMD-structured coordinates format, and forceTorque
namespace uammd{
namespace structured{

//GROMACS convention
//rij = rj - ri
//Fij force over particle i due to j
//Virial -real(0.5)*dot(rij,Fij)
//Stress -real(0.5)*outer(rij,Fij)

  __device__ __host__ real computeVirial(const real3& rij,const real3& Fij){
    return real(0.5)*dot(-rij,Fij);
  }

  __device__ __host__ tensor3 computeStress(const real3& rij,const real3& Fij){
    return real(0.5)*outer(-rij,Fij);
  }

  __device__ __host__ tensor3 computeHessianRadialPotential(const real3& rij, const real& invr,
							    const real& invr2, const real& dudr,
							    const real& d2udr2){

    tensor3 H = tensor3(0.0);

    const real rxx = rij.x*rij.x*invr2;
    const real rxy = rij.x*rij.y*invr2;
    const real rxz = rij.x*rij.z*invr2;
    const real ryy = rij.y*rij.y*invr2;
    const real ryz = rij.y*rij.z*invr2;
    const real rzz = rij.z*rij.z*invr2;

    H.xx = invr*(real(-1.0)+rxx)*dudr-rxx*d2udr2;
    H.xy = rxy*invr*dudr - rxy*d2udr2;
    H.xz = rxz*invr*dudr - rxz*d2udr2;

    H.yx = H.xy;
    H.yy = invr*(real(-1.0)+ryy)*dudr-ryy*d2udr2;
    H.yz = ryz*invr*dudr - ryz*d2udr2;

    H.zx = H.xz;
    H.zy = H.yz;
    H.zz = invr*(real(-1.0)+rzz)*dudr-rzz*d2udr2;
    return H;
  }

  __device__ __host__ tensor3 computeGradientAngle2AngularPotential(const real3& rji, const real3& rjk,
								    const real3& rik, const real& invrji,
								    const real& invrjk, const real& invrik,
								    const real& invrji2, const real& invrjk2,
								    const real& invrik2, const real3& gradTheta1,
								    const real3& gradTheta2,
								    const real& sijk, const real& cijk,
								    const int& id1, const int& id2){

    tensor3 grad2theta;

    tensor3 I = tensor3();
    I.xx = real(1.0);
    I.yy = real(1.0);
    I.zz = real(1.0);

    tensor3 gradtheta2_gradtheta1 = outer(gradTheta2,gradTheta1);
    real cotanijk = cijk/sijk;

    if (id1 == 0 and id2 == 0){ //ii
      tensor3 term1 = cotanijk*(-real(1.0)*outer(rji,rji)*invrji2*invrji2-gradtheta2_gradtheta1+I*invrji2);
      tensor3 term2 = -(outer(gradTheta1,rji)+outer(rji,gradTheta1))*invrji2;
      grad2theta = term1 + term2;

    } else if (id1 == 0 and id2 == 1){//ij, ji
      tensor3 term1 = outer(-cotanijk*gradTheta2 + invrji2*rji + invrjk2*rjk, gradTheta1);
      tensor3 term2 = I*(1-invrji/invrjk*cijk);
      tensor3 term3 = cijk*invrji*invrjk*outer(invrji2/invrjk2*rji-rjk,rji);
      tensor3 term4 = -invrji/invrjk*sijk*outer(gradTheta2,rji);

      tensor3 kk = outer(rjk, gradTheta2);

      grad2theta = term1+invrji*invrjk/sijk*(term2+term3+term4);

    } else if (id1 == 0 and id2 == 2){//ik, ki
      tensor3 term1 = cotanijk*(-real(1.0)*gradtheta2_gradtheta1+invrji2*invrjk2*outer(rjk,rji));
      tensor3 term2 = -(invrjk2*outer(rjk,gradTheta1)+invrji2*outer(gradTheta2,rji))-invrjk*invrji/sijk*I;
      grad2theta = term1 + term2;


    } else if (id1 == 1 and id2 == 1){//jj
      tensor3 term1 = -cotanijk*gradtheta2_gradtheta1;
      tensor3 term2 = (-2*(invrji*invrjk)+cijk*(invrji2+invrjk2))*I;
      tensor3 term3 = outer(rjk*invrjk2+rji*invrji2, (rji+rjk)*invrji*invrjk);
      tensor3 term4 = -2*cijk*(outer(rji,rji)*invrji2*invrji2+outer(rjk,rjk)*invrjk2*invrjk2);
      tensor3 term5 = sijk*outer(gradTheta2, rjk*invrjk2+rji*invrji2);
      grad2theta = term1 + (term2 + term3 + term4 + term5)/sijk;

    } else if (id1 == 1 and id2 == 2){//jk, kj
      tensor3 term1 = outer(gradTheta2, -cotanijk*gradTheta1 + invrji2*rji + invrjk2*rjk);
      tensor3 term2 = I*(1-invrjk/invrji*cijk);
      tensor3 term3 = cijk*invrjk*invrji*outer(rjk, invrjk2/invrji2*rjk-rji);
      tensor3 term4 = -invrjk/invrji*sijk*outer(rjk, gradTheta1);

      tensor3 kk = outer(rjk, gradTheta1);
      grad2theta = term1+invrji*invrjk/sijk*(term2+term3+term4);

    } else if (id1 == 2 and id2 == 2){//kk
      tensor3 term1 = cotanijk*(-real(1.0)*outer(rjk,rjk)*invrjk2*invrjk2-gradtheta2_gradtheta1+I*invrjk2);
      tensor3 term2 = -(outer(gradTheta1,rjk)+outer(rjk,gradTheta1))*invrjk2;
      grad2theta = term1 + term2;

    }
    return grad2theta;
  }

SFINAE_DEFINE_HAS_MEMBER(getEnergyTransverser);
SFINAE_DEFINE_HAS_MEMBER(getForceTransverser);
SFINAE_DEFINE_HAS_MEMBER(getLambdaTransverser);
SFINAE_DEFINE_HAS_MEMBER(getStressTransverser);
SFINAE_DEFINE_HAS_MEMBER(getForceTorqueMagneticFieldTransverser);
SFINAE_DEFINE_HAS_MEMBER(getMagneticFieldTransverser);
SFINAE_DEFINE_HAS_MEMBER(getHessianTransverser);
SFINAE_DEFINE_HAS_MEMBER(getPairwiseForceTransverser);

}}

///////////////////////////////////////////////

//Groups List (needed by Input)
#include"DataStructures/GroupsList/GroupsList.cuh"

//Input from file
#include"InputOutput/Input/InputFormats/InputJSON.cuh"

//System
#include"System/ExtendedSystem.cuh"

//Input manager
#include"InputOutput/Input/Input.cuh"

//GlobalData
#include"GlobalData/GlobalData.cuh"

//ParticleData
#include"ParticleData/StateLoader.cuh"
#include"ParticleData/ExtendedParticleData.cuh"

//GroupsUtils
#include"ParticleGroup/ParticleGroupUtils.cuh"

//NeighbourList
#include"DataStructures/ExclusionsList/ExclusionsList.cuh"
#include"DataStructures/VerletConditionalListSet/VerletConditionalListSetBase.cuh"
#include"DataStructures/VerletConditionalListSet/VerletConditionalListSetFactory.cuh"
#include"DataStructures/VerletConditionalListSet/VerletConditionalListSetLoaders.cuh"
#include"DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cuh"

//Parameter handler
#include"Utils/ParameterHandler/CheckDataConsistency.cuh"
#include"Utils/ParameterHandler/SingleParameterHandler.cuh"
#include"Utils/ParameterHandler/PairParameterHandler.cuh"

//Output formats
#include"InputOutput/Output/Output.cuh"

//Measures
#include"Utils/Measures/MeasuresBasic.cuh"

//Interactors
#include"Interactor/Single/SingleInteractor.cuh"
#include"Interactor/Pair/PairInteractor.cuh"
#include"Interactor/Bonds/BondsInteractor.cuh"
#include"Interactor/Set/SetInteractor.cuh"
#include"Interactor/AFM/AFMInteractor.cuh"

//Potentials
#include"Interactor/BasicPotentials.cuh"
#include"Interactor/BasicParameters.cuh"

//PatchyParticles
#include"Interactor/Patches/PatchesIncluders.cuh"
#include"Interactor/Patches/GenericPatchesPotentialLoader.cuh"
#include"Interactor/PatchyParticles/PatchyParticlesInteractor.cuh"
#include"Interactor/PatchyParticles/PatchyParticlesIncluders.cuh"

#include"Interactor/Bonds/BondsIncluders.cuh"
#include"Interactor/Single/SingleIncluders.cuh"
#include"Interactor/Set/SetIncluders.cuh"
#include"Interactor/Pair/PairIncluders.cuh"
#include"Interactor/AFM/AFMIncluders.cuh"

#include"Interactor/GenericPotentialLoader.cuh"

//Topology
#include"Topology/Topology.cuh"

//Force fields
#include"ForceFields/ForceFields.cuh"

//Backup
#include"Utils/Backup/Backup.cuh"

//Integrators
#include"Integrator/IntegratorBase.cuh"
//#include"Integrator/Minimization/MinimizationIncluders.cuh"
//#include"Integrator/Special/SpecialIncluders.cuh"
//#include"Integrator/NVE/NVEIncluders.cuh"
//#include"Integrator/NVT/NVTIncluders.cuh"
#include"Integrator/IntegratorFactory.cuh"
#include"Integrator/IntegratorLoaders.cuh"
#include"Integrator/IntegratorManager.cuh"

//Simulations
#include"SimulationStep/SimulationStep.cuh"
#include"SimulationStep/SimulationStepFactory.cuh"

#include"Utils/Backup/BackupStep.cuh"

#include"SimulationStep/GenericSimulationStepLoader.cuh"
#include"SimulationStep/SimulationStepManager.cuh"

#include"Simulation/Simulation.cuh"

#endif
