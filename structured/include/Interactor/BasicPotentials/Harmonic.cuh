#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    struct Harmonic{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                const real& K,const real& r0){

            const real r        = sqrt(r2);
            const real r0_div_r = ((r0 == r)?real(1.0):r0/r);
            const real fmod     = K*(real(1.0)-r0_div_r);

            return fmod*rij;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                const real& K,const real& r0){

            const real dr = sqrt(r2)-r0;

            return real(0.5)*K*dr*dr;

        }

      static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
                                               const real& K,const real& r0){
	    
        tensor3 H = tensor3(0.0);
        
        const real invr2 = real(1.0)/r2;
        const real invr  = sqrt(invr2);
	    
        H.xx =  (r0 == real(0.0))?-K:K*(real(-1.0)+r0*invr*(real(1.0)-rij.x*rij.x*invr2));
        H.xy =  (r0 == real(0.0))?real(0.0):-K*r0*invr*rij.x*rij.y*invr2;
        H.xz =  (r0 == real(0.0))?real(0.0):-K*r0*invr*rij.x*rij.z*invr2;
        
        H.yx =  H.xy;
        H.yy =  (r0 == real(0.0))?-K:K*(real(-1.0)+r0*invr*(real(1.0)-rij.y*rij.y*invr2));
        H.yz =  (r0 == real(0.0))?real(0.0):-K*r0*invr*rij.y*rij.z*invr2;
        
        H.zx = H.xz;
        H.zy = H.yz;
        H.zz = (r0 == real(0.0))?-K:K*(real(-1.0)+r0*invr*(real(1.0)-rij.z*rij.z*invr2));
        
        return H;
      }
    };

    struct lambdaHarmonic{

        static inline __device__ real3 force(const real3& rij,   const real& r2,
                const real& K,      const real& r0,
                const real& lambda, const int& n){

            const real lambda_n = pow(lambda,n);

            const real r        = sqrt(r2);
            const real r0_div_r = ((r0 == r)?real(1.0):r0/r);

            const real fmod = lambda_n*K*(real(1.0)-r0_div_r);

            return fmod*rij;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                const real& K,    const real& r0,
                const real& lambda, const int& n){

            const real dr       = sqrt(r2)-r0;
            const real lambda_n = pow(lambda,n);

            return real(0.5)*lambda_n*K*dr*dr;

        }

        static inline __device__ real lambdaDerivative(const real3& rij, const real& r2,
                const real& K,    const real& r0,
                const real& lambda, const int& n){
            const real dr = sqrt(r2)-r0;

            const real lambda_n_der = n*pow(lambda,n-1);

            return real(0.5)*lambda_n_der*K*dr*dr;

        }
    };

    struct HarmonicAnisotropic{

        static inline __device__ real3 force(const real3& rij,const real3& K,const real3& r0){
            real3 rij_abs = abs(rij);
            real3 sign_rij;
            sign_rij.x = ((rij.x > real(0.0))?real(1.0):real(-1.0));
            sign_rij.y = ((rij.y > real(0.0))?real(1.0):real(-1.0));
            sign_rij.z = ((rij.z > real(0.0))?real(1.0):real(-1.0));
            return {K.x*(rij_abs.x-r0.x)*sign_rij.x,
                K.y*(rij_abs.y-r0.y)*sign_rij.y,
                K.z*(rij_abs.z-r0.z)*sign_rij.z};
        }

        static inline __device__ real energy(const real3& rij,const real3& K,const real3& r0){
            real3 rij_abs = abs(rij);
            return real(0.5)*K.x*(rij_abs.x-r0.x)*(rij_abs.x-r0.x)+
                real(0.5)*K.y*(rij_abs.y-r0.y)*(rij_abs.y-r0.y)+
                real(0.5)*K.z*(rij_abs.z-r0.z)*(rij_abs.z-r0.z);

        }

    };

    struct lambdaHarmonicAnisotropic{

        static inline __device__ real3 force(const real3& rij,
                const real3& K,const real3& r0,
                const real& lambda, const int& n){

            const real lambda_n = pow(lambda,n);
            real3 rij_abs       = abs(rij);
            real3 sign_rij;
            sign_rij.x = ((rij.x > real(0.0))?real(1.0):real(-1.0));
            sign_rij.y = ((rij.y > real(0.0))?real(1.0):real(-1.0));
            sign_rij.z = ((rij.z > real(0.0))?real(1.0):real(-1.0));

            return {lambda_n*K.x*(rij_abs.x-r0.x)*sign_rij.x,
                lambda_n*K.y*(rij_abs.y-r0.y)*sign_rij.y,
                lambda_n*K.z*(rij_abs.z-r0.z)*sign_rij.z};
        }

        static inline __device__ real energy(const real3& rij,
                const real3& K,const real3& r0,
                const real& lambda, const int& n){

            const real lambda_n = pow(lambda,n);
            real3 rij_abs       = abs(rij);
            return real(0.5)*lambda_n*K.x*(rij_abs.x-r0.x)*(rij_abs.x-r0.x)+
                real(0.5)*lambda_n*K.y*(rij_abs.y-r0.y)*(rij_abs.y-r0.y)+
                real(0.5)*lambda_n*K.z*(rij_abs.z-r0.z)*(rij_abs.z-r0.z);

        }

        static inline __device__ real lambdaDerivative(const real3& rij,
                const real3& K,const real3& r0,
                const real& lambda, const int& n){
            const real lambda_n_der = n*pow(lambda,n-1);
            real3 rij_abs       = abs(rij);

            return real(0.5)*lambda_n_der*K.x*(rij_abs.x-r0.x)*(rij_abs.x-r0.x)+
                real(0.5)*lambda_n_der*K.y*(rij_abs.y-r0.y)*(rij_abs.y-r0.y)+
                real(0.5)*lambda_n_der*K.z*(rij_abs.z-r0.z)*(rij_abs.z-r0.z);
        }

    };

}}}}
