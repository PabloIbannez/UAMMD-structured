#pragma once

#include "utils/container.h"
#include "utils/quaternion.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{
namespace Magnetic{

    template<class T>  using cached_vector = uninitialized_cached_vector<T>;

  namespace LLG{
    /* Computes dm/dt = -gyroRation/(1+damping^2)*(m x b+ damping * m x m x b)*/
    inline __device__ real3 computeMagnetizationDerivative(real3 field, real3 magnetization,
							   real damping, real gyroRatio){
	real3 mxbi = cross(magnetization,field);
	real3 amxmxbi = damping*cross(magnetization,mxbi);
	real3 dmi = -gyroRatio/(1+damping*damping)*(mxbi+amxmxbi);
	return dmi;
      }

      /* Computes the thermal field i.e. an hypothetial field to consider the thermal fluctuations
       * in the system. The noise follow a Gaussian distribution with mean 0 and standard deviation
       * dW = sqrt(2*k_BT*damping/gyroRation*msat*V*dt)
       */
      inline __device__ real3 computeThermalField(real fluctuationsAmplitude,
						  int step, uint seed, int id){
	Saru rng(seed,id,step);
	real3 randomField = make_real3(rng.gf(0, fluctuationsAmplitude),
				       rng.gf(0, fluctuationsAmplitude).x);
	return randomField;
      }

      /* Computes the anisotropy field. i.e. an internal field that pushes the magnetization
       * towards the direction of the anisotropy axis. B_anis = 2*anisCons/msat*(m·u)u
       */
      inline __device__ real3 computeAnisotropyField(real3 mi, real anisConst_div_msat){
	real dotm_u = mi.z; //In the particles frame the easy axis is in the z axis
	real fieldAnisParticle_z = real(2.0)*anisConst_div_msat*dotm_u;
	return {0, 0, fieldAnisParticle_z};
      }

      inline __device__ real3 computeAnisotropyField(Quat diri, real3 mi, real anisConst_div_msat){
	real3 vz = diri.getVz();
	real dotm_u = dot(mi, vz); //In the particles frame the easy axis is in the z axis
	real3 fieldAnisParticle = real(2.0)*anisConst_div_msat*dotm_u*vz;
	return fieldAnisParticle;
      }

      namespace Euler{

	inline __global__ void integrateLLG(real4* dir, real4*  field, real* anisotropy,
				     real4* magnetization,
				     ParticleGroup::IndexIterator indexIterator,
				     real dt, real kbT, real damping, real msat, real gyroRatio,
				     int currentStep, int seed,  int N){

	  int id = blockIdx.x*blockDim.x+threadIdx.x;
	  if(id>=N) return;
	  int i = indexIterator[id];
	  real4 m_and_M = magnetization[i];
	  real3 mi = make_real3(m_and_M);
	  real Mi = m_and_M.w;
	  real anisotropy_i = anisotropy[i];
	  if (Mi == real(0.0)) return;
	  real3 bi = make_real3(field[i]);
	  //Moves from the frame of the laboratory to the frame of the particle.
	  bi = rotateVector(Quat(dir[i]).getConjugate(), bi);
	  real fluctuationsAmplitude = sqrt(2*kbT*damping/(gyroRatio*Mi*dt));
	  bi+=computeThermalField(fluctuationsAmplitude, currentStep, seed, id);
	  real3 anisotropyField = computeAnisotropyField(mi, anisotropy_i/msat);
	  bi+=anisotropyField;
	  real3 dmi = computeMagnetizationDerivative(bi, mi, damping, gyroRatio)*dt;
	  mi += dmi;
	  magnetization[i] = make_real4(mi*rsqrt(dot(mi,mi)), Mi);
	}

	inline void updateMagnetization(std::shared_ptr<ParticleData> pd,
				 std::shared_ptr<ParticleGroup> pg,
				 std::shared_ptr<GlobalData> gd,
				 cudaStream_t st, real damping,
				 real msat, real gyroRatio,
				 real dt, real kBT){

	  auto field = pd->getMagneticField(access::location::gpu, access::mode::read).raw();
	  auto dir = pd->getDir(access::location::gpu, access::mode::read).raw();
	  auto anisotropy = pd->getAnisotropy(access::location::gpu, access::mode::read).raw();
	  auto magnetization = pd->getMagnetization(access::location::gpu, access::mode::readwrite).raw();
	  auto groupIterator = pg->getIndexIterator(access::location::gpu);
	  uint currentStep = gd->getFundamental()->getCurrentStep();
	  uint seed = gd->getSystem()->getSeed();
	  int BLOCKSIZE = 128;
	  int numberParticles = pg->getNumberParticles();
	  uint Nthreads = BLOCKSIZE<numberParticles?BLOCKSIZE:numberParticles;
	  uint Nblocks = numberParticles/Nthreads +  ((numberParticles%Nthreads!=0)?1:0);
	  integrateLLG<<<Nblocks, Nthreads, 0, st>>>(dir, field, anisotropy,
						     magnetization, groupIterator, dt,
						     kBT, damping, msat, gyroRatio,
						     currentStep, seed, numberParticles);
	}
      }

      namespace Heun{

	inline __global__ void integrateLLG(real4* dir, real4* field, real4* magnetization,
				     real3* initialMagnetization, real* anisotropy,
				     ParticleGroup::IndexIterator indexIterator,
				     real dt, real kbT, real damping, real msat,
				     real gyroRatio, int currentStep, int seed, int N,
				     int step){

	  int id = blockIdx.x*blockDim.x+threadIdx.x;
	  if(id>=N) return;
	  int i = indexIterator[id];
	  real4 m_and_M = magnetization[i];
	  real3 mi = make_real3(m_and_M);
	  real Mi = m_and_M.w;
	  if (Mi == real(0.0)) return;
	  Quat diri = dir[i];
	  real anisotropy_i = anisotropy[i];
	  if (step==0){
	    initialMagnetization[i] = mi;
	  }
	  real3 bi = make_real3(field[i]);
	  //Moves from the frame of the laboratory to the frame of the particle.
	  bi = rotateVector(Quat(dir[i]).getConjugate(), bi);
	  real fluctuationsAmplitude = sqrt(2*kbT*damping/(gyroRatio*Mi*dt));
	  bi += computeThermalField(fluctuationsAmplitude, currentStep, seed, id);
	  real3 anisotropyField = computeAnisotropyField(mi, anisotropy_i/msat);
	  bi+=anisotropyField;
	  real3 dmi = computeMagnetizationDerivative(bi,mi,damping, gyroRatio)*dt;
	  if (step==0){
	    mi+=dmi;
	    field[i] = real4();
	  } else if (step==1) {
	    real3 mi_t0 = initialMagnetization[i];
	    mi = real(0.5)*(mi_t0 + mi + dmi);
	    mi *= rsqrt(dot(mi,mi));
	  }
	  magnetization[i] = make_real4(mi, Mi);
	}

	inline void updateHalfStep(std::shared_ptr<ParticleData> pd,
			    std::shared_ptr<ParticleGroup> pg,
			    std::shared_ptr<GlobalData> gd,
			    cached_vector<real3> &initialMagnetization,
			    real damping, real msat, real gyroRatio,
			    real dt, real kBT, int step, cudaStream_t st){


	  auto field = pd->getMagneticField(access::location::gpu, access::mode::read).raw();
	  auto dir = pd->getDir(access::location::gpu, access::mode::read).raw();
	  auto anisotropy = pd->getAnisotropy(access::location::gpu, access::mode::read).raw();
	  auto magnetization = pd->getMagnetization(access::location::gpu, access::mode::readwrite).raw();
	  auto groupIterator = pg->getIndexIterator(access::location::gpu);
	  auto initMagnet_ptr = thrust::raw_pointer_cast(initialMagnetization.data());
	  uint currentStep = gd->getFundamental()->getCurrentStep();
	  uint seed = gd->getSystem()->getSeed();
	  int BLOCKSIZE = 128;
	  int numberParticles = pg->getNumberParticles();
	  uint Nthreads = BLOCKSIZE<numberParticles?BLOCKSIZE:numberParticles;
	  uint Nblocks = numberParticles/Nthreads +  ((numberParticles%Nthreads!=0)?1:0);
	  integrateLLG<<<Nblocks, Nthreads, 0, st>>>(dir, field, magnetization,
						     initMagnet_ptr, anisotropy,
						     groupIterator, dt, kBT,
						     damping, msat, gyroRatio,
						     currentStep, seed,
						     numberParticles, step);
	}
      }
    }

}}}}}
