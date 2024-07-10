#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct Contact{

        //Force
        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            real eps;
            if(epsilon==real(0.0)){
                eps=  fabs(zeroEnergy);
            } else {
                eps = fabs(epsilon);
            }

            real3 fs = Steric::Steric::force<12>(rij,r2,eps,sigma);

            real r     = sqrt(r2);
            real sech2 = real(1.0)/cosh(n*(r0-r));
                 sech2 = sech2*sech2;

            real3 f  = -real(0.5)*epsilon*n*sech2*rij/r;

            return fs+f;
        }

        //Virial
        static inline __device__ real virial(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            return real(0);
        }

        //Stress
        static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                const real& epsilon,const real& sigma,
                                                const real& n,const real& r0,
                                                const real& zeroEnergy){

            return tensor3(0);
        }

        //Energy
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            real eps;
            if(epsilon==real(0.0)){
                eps=zeroEnergy;
            } else {
                eps = epsilon;
            }

            real es = Steric::Steric::energy<12>(rij,r2,eps,sigma);
            real e  = real(0.5)*epsilon*(real(1.0)+tanh(n*(r0-sqrt(r2))));

            return es+e;
        }
    };

}}}}
