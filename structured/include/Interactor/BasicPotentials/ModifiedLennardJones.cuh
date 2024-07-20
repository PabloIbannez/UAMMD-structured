#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/WCA.cuh"
#include "Interactor/BasicPotentials/GaussianWell.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    namespace ModifiedLennardJones{

        struct Gaussian{

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){

                return WCA::Type2::force(rij,r2,epsilon,sigma)+
                       GaussianWell::force(rij,r2,epsilon,sigma,D);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma,const real& D){
                return tensor3(0);
            }

	  //Energy
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						     const real& epsilon,const real& sigma,const real& D){
                return WCA::Type2::hessian(rij,r2,epsilon,sigma)+
                       GaussianWell::hessian(rij,r2,epsilon,sigma,D);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return WCA::Type2::energy(rij,r2,epsilon,sigma)+
                       GaussianWell::energy(rij,r2,epsilon,sigma,D);
            }

        };
    }

}}}}
