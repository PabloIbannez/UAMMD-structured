#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct GaussianWell{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                const real& e,const real& r0,const real& D){

            const real r = sqrt(r2);
            const real dr = r-r0;

            const real fmod = (e/D)*(real(1.0)-r0/r)*exp(-dr*dr/(real(2.0)*D));
            return fmod*rij;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                const real& e,const real& r0,const real& D){

            const real dr = sqrt(r2)-r0;

            return -e*exp(-dr*dr/(real(2.0)*D))/real(2.0);

        }

        static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
                const real& e,const real& r0,const real& D){

            const real r  = sqrt(r2);
            const real dr = r-r0;

            const real dudr   = (e/D)*(real(1.0)-r0/r)*exp(-dr*dr/(real(2.0)*D));
            const real d2udr2 = real(0.5)*e*(D - dr*dr)*exp(-real(0.5)*dr*dr/D)/(D*D);

            return computeHessianRadialPotential(rij, real(1.0)/r, real(1.0)/r2, dudr, d2udr2);

        }

    };

}}}}
