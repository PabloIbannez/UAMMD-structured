#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/StiffnessMask.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct DistanceSwitchCosine{

        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& rc,const real& K){

            const real r = sqrtf(r2);

            real  swt;

            if(r > rc){
                swt = real(1.0);
            } else {
                swt = cosf(real(M_PI)*r/rc);
                swt = real(0.5)*(real(1.0)+swt);
                swt = StiffnessMask::stiffnessMask(swt,K);
            }

            return swt;
        }

        static inline __device__ real3 force(const real3& rij,const real& r2,
                                             const real& rc  ,const real& K){

            const real r = sqrtf(r2);

            if(r < real(1e-6)){
                return make_real3(0.0);
            }

            real  swt;
            real  swt_der;

            if(r > rc){
                swt_der = real(0.0);
            } else {
                swt = cosf(real(M_PI)*r/rc);
                swt = real(0.5)*(real(1.0)+swt);

                swt_der = -(real(0.5)*real(M_PI)/rc)*sinf(real(M_PI)*r/rc);

                swt_der = StiffnessMask::stiffnessMaskFirstDerivative(swt,K)*swt_der;
            }

            const real3 f = swt_der*(rij/r);

            return f;
        }

        static inline __device__ real4 forceEnergy(const real3& rij, const real& r2,
                                                   const real& rc,const real& K){
            const real r = sqrtf(r2);

            if(r < real(1e-6)){
                return make_real4(0.0);
            }

            real  e = energy(rij,r2,rc,K);
            real3 f = force(rij,r2,rc,K);

            return make_real4(f,e);
        }

    };

}}}}
