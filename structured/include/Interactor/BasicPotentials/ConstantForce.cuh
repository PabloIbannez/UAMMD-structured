#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct ConstantForce{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& F){
            const real invr = rsqrt(r2);
            return -F*rij*invr;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& F){
            return sqrtf(r2)*F;

        }
    };

}}}}
