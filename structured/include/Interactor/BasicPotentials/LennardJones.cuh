#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/Steric.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    namespace LennardJones{

        //U(r)=4e((s/r)^12-(s/r)^6)
        struct Type1{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(24.0)*epsilon*(real(2.0)*sinvr12-sinvr6)*invr2;

                return fmod*rij;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real e = real(4.0)*epsilon*(sinvr12-sinvr6);

                return e;
            }


	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon, const real& sigma){

	    const real  invr2  = real(1.0)/r2;
	    const real  invr   = sqrt(invr2);
	    const real sinvr2  = sigma*sigma*invr2;
	    const real sinvr6  = sinvr2*sinvr2*sinvr2;
	    const real sinvr12 = sinvr6*sinvr6;

	    const real   dudr =  real(24.0)*epsilon*(sinvr6 - real(2.0)*sinvr12)*invr;
	    const real d2udr2 = -real(24.0)*epsilon*(real(7.0)*sinvr6 - real(26.0)*sinvr12)*invr2;

	    return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
        }

        };

        //U(r)=e((s/r)^12-2*(s/r)^6)
        struct Type2{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(12.0)*epsilon*(sinvr12-sinvr6)*invr2;

                return fmod*rij;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(sinvr12-real(2.0)*sinvr6);

                return e;
            }

	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon, const real& sigma){

	    const real  invr2  = real(1.0)/r2;
	    const real  invr   = sqrt(invr2);
	    const real sinvr2  = sigma*sigma*invr2;
	    const real sinvr6  = sinvr2*sinvr2*sinvr2;
	    const real sinvr12 = sinvr6*sinvr6;

	    const real   dudr =  real(12.0)*epsilon*(sinvr6 - sinvr12)*invr;
	    const real d2udr2 = -real(12.0)*epsilon*(real(7.0)*sinvr6 - real(13.0)*sinvr12)*invr2;

	    return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
	  }

        };

        //U(r)=e(5(s/r)^12-6(s/r)^10)
        struct Type3{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(60.0)*epsilon*(sinvr12-sinvr10)*invr2;

                return fmod*rij;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(real(5.0)*sinvr12-real(6.0)*sinvr10);

                return e;
            }

	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon, const real& sigma){

	    const real  invr2  = real(1.0)/r2;
	    const real  invr   = sqrt(invr2);
	    const real sinvr2  = sigma*sigma*invr2;
	    const real sinvr4  = sinvr2*sinvr2;
	    const real sinvr10 = sinvr4*sinvr4*sinvr2;
	    const real sinvr12 = sinvr10*sinvr2;

	    const real   dudr =  real(60.0)*epsilon*(sinvr10 - sinvr12)*invr;
	    const real d2udr2 = -real(60.0)*epsilon*(real(11.0)*sinvr10 - real(13.0)*sinvr12)*invr2;

	    return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
	  }

        };

        //U(r)=e(13(s/r)^12-18(s/r)^10+4(s/r)^6)
        struct KaranicolasBrooks{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -epsilon*(real(156.0)*sinvr12-real(180.0)*sinvr10+real(24.0)*sinvr6)*invr2;

                return fmod*rij;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(real(13.0)*sinvr12-real(18.0)*sinvr10+real(4.0)*sinvr6);

                return e;
            }

        };

        struct SoftCoreType1{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                //I assume compiler will optimize this
                const real3 f = Steric::Steric12SoftCore::force(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n)
                               -Steric::Steric6SoftCore::force(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n);

                return f;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                return computeVirial(rij,force(rij,r2,epsilon,sigma,lambda,alpha,n));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma,
                                                    const real& lambda, const real& alpha,
                                                    const int& n){

                return computeStress(rij,force(rij,r2,epsilon,sigma,lambda,alpha,n));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                const real e = Steric::Steric12SoftCore::energy(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n)
                              -Steric::Steric6SoftCore::energy(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n);

                return e;
            }

            static inline __device__ real lambdaDerivative(const real3& rij, const real& r2,
                                                           const real& epsilon,const real& sigma,
                                                           const real& lambda, const real& alpha,
                                                           const int& n){

                const real ld = Steric::Steric12SoftCore::lambdaDerivative(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n)
                               -Steric::Steric6SoftCore::lambdaDerivative(rij,r2,real(4.0)*epsilon,sigma,lambda,alpha,n);

                return ld;
            }

        };

        struct SoftCoreType2{

            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                //I assume compiler will optimize this
                const real3 f = Steric::Steric12SoftCore::force(rij,r2,epsilon,sigma,lambda,alpha,n)
                               -real(2.0)*Steric::Steric6SoftCore::force(rij,r2,epsilon,sigma,lambda,alpha,n);

                return f;
            }

            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                return computeVirial(rij,force(rij,r2,epsilon,sigma,lambda,alpha,n));
            }

            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma,
                                                    const real& lambda, const real& alpha,
                                                    const int& n){

                return computeStress(rij,force(rij,r2,epsilon,sigma,lambda,alpha,n));
            }

            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                const real e = Steric::Steric12SoftCore::energy(rij,r2,epsilon,sigma,lambda,alpha,n)
                              -real(2.0)*Steric::Steric6SoftCore::energy(rij,r2,epsilon,sigma,lambda,alpha,n);

                return e;
            }

            static inline __device__ real lambdaDerivative(const real3& rij, const real& r2,
                                                           const real& epsilon,const real& sigma,
                                                           const real& lambda, const real& alpha,
                                                           const int& n){

                const real ld = Steric::Steric12SoftCore::lambdaDerivative(rij,r2,epsilon,sigma,lambda,alpha,n)
                               -real(2.0)*Steric::Steric6SoftCore::lambdaDerivative(rij,r2,epsilon,sigma,lambda,alpha,n);

                return ld;
            }

        };

    }

}}}}
