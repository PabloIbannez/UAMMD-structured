#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct Fene{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                const real& K,const real& r0,const real& R0){

            const real R02 = R0*R0;

            const real  r  = sqrt(r2);
            const real dr  = r-r0;
            const real dr2 = dr*dr;

            const real fmod = K*dr/(real(1.0)-dr2/R02);

            return fmod*rij/r;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                const real& K,const real& r0,const real& R0){

            const real R02 = R0*R0;

            const real  r  = sqrt(r2);
            const real dr  = r-r0;
            const real dr2 = dr*dr;

            return -(K*R02/real(2.0))*log(real(1.0)-dr2/R02);

        }

        static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
                const real& K,const real& r0, const real& R0){
            const real R02    = R0*R0;
            const real r      = sqrt(r2);
            const real invr   = real(1.0)/r;
            const real invr2  = invr*invr;
            const real dr     = r-r0;
            const real dr2    = dr*dr;

            const real dudr   = K*dr/(real(1.0)-dr2/R02);
            const real d2udr2 = K*R02*(R02+dr2)/((R02-dr2)*(R02-dr2));

            return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
        }
    };

}}}}
