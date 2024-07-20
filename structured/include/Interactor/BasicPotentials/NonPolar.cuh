#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/Steric.cuh"
#include "Interactor/BasicPotentials/LennardJones.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{


    struct NonPolar{

        //Force
        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){

            if(epsilon == real(0.0) ){
                return BasicPotentials::Steric::Steric::force<12>(rij,r2,zeroEnergy,sigma);
            }

            real3 f = BasicPotentials::LennardJones::Type1::force(rij,r2,fabs(epsilon),sigma);
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0) and r2>=r02){
                f=-f;
            }

            return f;
        }

        //Virial
        static inline __device__ real virial(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
            return real(0);
        }

        //Stress
        static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                const real& epsilon,const real& sigma,const real& zeroEnergy){
            return tensor3(0);
        }

      //Hessian
      static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
					       const real& epsilon,const real& sigma,const real& zeroEnergy){
	if(epsilon == real(0.0) ){
	  return BasicPotentials::Steric::Steric::hessian<12>(rij,r2,zeroEnergy,sigma);
	}

            tensor3 H = BasicPotentials::LennardJones::Type1::hessian(rij,r2,fabs(epsilon),sigma);
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0) and r2>=r02){
                H=-H;
            }

            return H;
      }

        //Energy
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){

            if(epsilon == real(0.0) ){
                return BasicPotentials::Steric::Steric::energy<12>(rij,r2,zeroEnergy,sigma);
            }

            real e  = BasicPotentials::LennardJones::Type1::energy(rij,r2,fabs(epsilon),sigma);

            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0)){

                if(r2<r02){
                    e+=real(2.0)*epsilon;
                } else {
                    e=-e;
                }
            }

            return e*real(0.5);
        }

    };

}}}}
