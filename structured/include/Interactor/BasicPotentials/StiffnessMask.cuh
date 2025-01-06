#pragma once

#include "uammd.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{
namespace StiffnessMask{

    // This StiffnessMask is tough to be used in combination with other function f(x) that:
    // 1. Return values in [0,1]
    // 2. Has a minimum at x=1
    //
    // Mainly StiffnessMask adds an additional parameter, K, that permits to
    // control the stiffness of the base function.

    static inline __device__ real stiffnessMask(const real& x,
                                                const real& K){
        return x/(x+K*(real(1.0)-x));
    }

    static inline __device__ real stiffnessMaskFirstDerivative(const real& x,
                                                               const real& K){
        const real denom = x+K*(real(1.0)-x);

        return K/(denom*denom);
    }

    static inline __device__ real stiffnessMaskSecondDerivative(const real& x,
                                                                const real& K){
        const real denom = x+K*(real(1.0)-x);

        return real(2.0)*K*(K-real(1.0))/(denom*denom*denom);
    }

    // Then F(x;K) = stiffnessMask(f(x),K)
    // F'(x;K) = stiffnessMaskFirstDerivative(f(x),K)*f'(x)
    // F''(x;K) = stiffnessMaskSecondDerivative(f(x),K)*f'(x)^2 + stiffnessMaskFirstDerivative(f(x),K)*f''(x)

}}}}}
