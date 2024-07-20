#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct DistanceSwitchExponential{

        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& E,const real& rc,const real& K){

            const real r = sqrtf(r2);
            const real A = expf(-real(2.0)*K);

            if(r >= rc){
                return real(0.0);
            }

            const real swt = expf(-K*(real(1.0)-cosf(real(M_PI)*r/rc)));

            return -E*(swt-A)/(real(1.0)-A);
        }

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& E,const real& rc,const real& K){
            const real r = sqrtf(r2);
            const real A = expf(-real(2.0)*K);

            if(r >= rc or r < real(1e-6)){
                return make_real3(0.0);
            }

            real swt = expf(-K*(real(1.0)-cosf(real(M_PI)*r/rc)));
                 swt = K*(real(M_PI)/rc)*sinf(real(M_PI)*r/rc)*swt;

            return E*(swt/(real(1.0)-A))*(rij/r);
        }

        static inline __device__ real4 forceEnergy(const real3& rij, const real& r2,
                                                   const real& E,const real& rc,const real& K){
            const real r = sqrtf(r2);
            const real A = expf(-real(2.0)*K);

            if(r >= rc or r < real(1e-6)){
                return make_real4(0.0,0.0,0.0,-E);
            }

            const real swt_e = expf(-K*(real(1.0)-cosf(real(M_PI)*r/rc)));
            const real swt_f = K*(real(M_PI)/rc)*sinf(real(M_PI)*r/rc)*swt_e;

            return make_real4(E*(swt_f/(real(1.0)-A))*(rij/r),-E*(swt_e-A)/(real(1.0)-A));
        }

    };

    struct DistanceSwitchCosine{

        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& E,real r_switch_start,real r_switch_end){

            const real r = sqrtf(r2);

            real  swt;

            if(r >= r_switch_end){
                swt = real(1.0);
            } else if(r < r_switch_start){
                swt = real(-1.0);
            } else {
                const real r_norm = (r-r_switch_start)/(r_switch_end-r_switch_start);
                swt = cosf(real(M_PI)*r_norm);

                swt = fminf(swt, real(1.0));
                swt = fmaxf(swt, real(-1.0));

                swt = -swt;

            }

            const real e = -E*(real(1.0)-swt)*real(0.5);

            return e;
        }

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& E,real r_switch_start,real r_switch_end){
            const real r = sqrtf(r2);

            if(r < real(1e-6)){
                return make_real3(0.0);
            }

            real  swt_fmod;

            if(r >= r_switch_end){
                swt_fmod = real(0.0);
            } else if(r < r_switch_start){
                swt_fmod = real(0.0);
            } else {
                const real r_norm = (r-r_switch_start)/(r_switch_end-r_switch_start);

                swt_fmod = sinf(real(M_PI)*r_norm);

                swt_fmod = fminf(swt_fmod, real(1.0));
                swt_fmod = fmaxf(swt_fmod, real(-1.0));

                swt_fmod = -real(M_PI)*swt_fmod/(r_switch_end-r_switch_start);
            }

            swt_fmod = -E*swt_fmod*real(0.5);

            const real3 f = swt_fmod*(rij/r);

            return f;
        }

        static inline __device__ real4 forceEnergy(const real3& rij, const real& r2,
                                                   const real& E,real r_switch_start,real r_switch_end){
            const real r = sqrtf(r2);

            if(r < real(1e-6)){
                return make_real4(0.0,0.0,0.0,-E);
            }

            real swt;
            real swt_fmod;

            if(r >= r_switch_end){
                swt = real(1.0);
                swt_fmod = real(0.0);
            } else if(r < r_switch_start){
                swt = real(-1.0);
                swt_fmod = real(0.0);
            } else {
                const real r_norm = (r-r_switch_start)/(r_switch_end-r_switch_start);

                swt = cosf(real(M_PI)*r_norm);

                swt = fminf(swt, real(1.0));
                swt = fmaxf(swt, real(-1.0));

                swt = -swt;

                swt_fmod = sinf(real(M_PI)*r_norm);

                swt_fmod = fminf(swt_fmod, real(1.0));
                swt_fmod = fmaxf(swt_fmod, real(-1.0));

                swt_fmod = -real(M_PI)*swt_fmod/(r_switch_end-r_switch_start);
            }

            swt      = -E*(real(1.0)-swt)*real(0.5);
            swt_fmod = -E*swt_fmod*real(0.5);

            const real3 f = swt_fmod*(rij/r);

            return make_real4(f,swt);
        }

    };

}}}}
