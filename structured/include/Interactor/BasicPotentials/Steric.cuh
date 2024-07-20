#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    namespace Steric{

        struct Steric{

            //Force
            template<int power>
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                          const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return make_real3(0);
            }

            //Virial
            template<int power>
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return real(0);
            }

            //Stress
            template<int power>
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return tensor3(0);
            }


            //Energy
            template<int power>
            static inline __device__ real energy(const real3& rij, const real& r2,
                                          const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return real(0);
            }

	  //Hessian
            template<int power>
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						     const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return tensor3(0);
            }


        };

        template<>
        inline __device__ real3 Steric::force<6>(const real3& rij, const real& r2,
                                         const real& epsilon,const real& sigma){

            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;

            const real fmod = -real(6.0)*epsilon*sinvr6*invr2;

            return fmod*rij;
        }

        template<>
        inline __device__ real3 Steric::force<12>(const real3& rij, const real& r2,
                                          const real& epsilon,const real& sigma){

            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;
            const real sinvr12 = sinvr6*sinvr6;

            const real fmod = -real(12.0)*epsilon*sinvr12*invr2;

            return fmod*rij;
        }


        template<>
        inline __device__ real Steric::virial<6>(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

            return computeVirial(rij,force<6>(rij,r2,epsilon,sigma));
        }

        template<>
        inline __device__ tensor3 Steric::stress<6>(const real3& rij, const real& r2,
                                            const real& epsilon,const real& sigma){

            return computeStress(rij,force<6>(rij,r2,epsilon,sigma));
        }

        template<>
        inline __device__ real Steric::virial<12>(const real3& rij, const real& r2,
                                                  const real& epsilon,const real& sigma){

            return computeVirial(rij,force<12>(rij,r2,epsilon,sigma));
        }

        template<>
        inline __device__ tensor3 Steric::stress<12>(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma){

            return computeStress(rij,force<12>(rij,r2,epsilon,sigma));
        }

        template<>
        inline __device__ real Steric::energy<6>(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;

            const real e = epsilon*sinvr6;

            return e;
        }

        template<>
        inline __device__ real Steric::energy<12>(const real3& rij, const real& r2,
                                          const real& epsilon,const real& sigma){

            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;
            const real sinvr12 = sinvr6*sinvr6;

            const real e = epsilon*sinvr12;

            return e;
        }

      template<>
      inline __device__ tensor3 Steric::hessian<6>(const real3& rij, const real& r2,
						   const real& epsilon, const real& sigma){

	const real  invr2  = real(1.0)/r2;
	const real   invr  = sqrt(invr2);
	const real sinvr2  = sigma*sigma*invr2;
	const real sinvr6  = sinvr2*sinvr2*sinvr2;

	const real dudr   = -real(6.0)*epsilon*sinvr6*invr;
	const real d2udr2 = real(-7.0)*dudr*invr;

	return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
      }

      template<>
      inline __device__ tensor3 Steric::hessian<12>(const real3& rij, const real& r2,
						   const real& epsilon, const real& sigma){

	const real  invr2   = real(1.0)/r2;
	const real   invr   = sqrt(invr2);
	const real sinvr2   = sigma*sigma*invr2;
	const real sinvr6   = sinvr2*sinvr2*sinvr2;
	const real sinvr12  = sinvr6*sinvr6;

	const real dudr   = real(-12.0)*epsilon*sinvr12*invr;
	const real d2udr2 = real(-13.0)*dudr*invr;

	return computeHessianRadialPotential(rij, invr, invr2, dudr, d2udr2);
      }

        struct Steric6{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                return Steric::force<6>(rij,r2,epsilon,sigma);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::virial<6>(rij,r2,epsilon,sigma);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return Steric::stress<6>(rij,r2,epsilon,sigma);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::energy<6>(rij,r2,epsilon,sigma);
            }

	   //Hessian
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						     const real& epsilon,const real& sigma){
                return Steric::hessian<6>(rij,r2,epsilon,sigma);
            }
        };

        struct Steric12{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                return Steric::force<12>(rij,r2,epsilon,sigma);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::virial<12>(rij,r2,epsilon,sigma);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return Steric::stress<12>(rij,r2,epsilon,sigma);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::energy<12>(rij,r2,epsilon,sigma);
            }

	        //Hessian
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return Steric::hessian<12>(rij,r2,epsilon,sigma);
            }

        };

        struct Steric6SoftCore{

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);

                const real denominator = alpha*lambdaFactor+rinvs6;

                return (lambda_n*epsilon/denominator);


            }

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);

                const real denominator = alpha*lambdaFactor+rinvs6;
                const real denominator2 = denominator*denominator;

                const real3 f = -epsilon*(lambda_n*real(6.0)/denominator2)*rinvs6*rij/r2;

                return f;
            }

            //Lambda derivative
            static inline __device__ real lambdaDerivative(const real3& rij, const real& r2,
                                                           const real& epsilon,const real& sigma,
                                                           const real& lambda, const real& alpha,
                                                           const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);
                const real lambda_n_der = n*pow(lambda,n-1);

                const real denominator = alpha*lambdaFactor+rinvs6;
                const real denominator2 = denominator*denominator;

                const real ld = epsilon*(lambda_n_der/denominator+(real(2.0)*lambda_n*alpha*(real(1.0)-lambda)/denominator2));

                return ld;
            }
        };

        struct Steric12SoftCore{

            //Energy
            static inline __device__ real energy(const real3& rij   , const real& r2,
                                                 const real& epsilon, const real& sigma,
                                                 const real& lambda , const real& alpha,
                                                 const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);

                const real denominator  = alpha*lambdaFactor+rinvs6;
                const real denominator2 = denominator*denominator;

                return (lambda_n*epsilon/denominator2);
            }

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,
                                                 const real& lambda, const real& alpha,
                                                 const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);

                const real denominator  = alpha*lambdaFactor+rinvs6;
                const real denominator2 = denominator*denominator;
                const real denominator3 = denominator2*denominator;

                const real3 f = -epsilon*(lambda_n*real(12.0)/denominator3)*rinvs6*rij/r2;

                return f;
            }

            //Lambda derivative
            static inline __device__ real lambdaDerivative(const real3& rij, const real& r2,
                                                           const real& epsilon,const real& sigma,
                                                           const real& lambda, const real& alpha,
                                                           const int& n){

                const real  invs2  = real(1.0)/(sigma*sigma);
                const real rinvs2  = r2*invs2;
                const real rinvs6  = rinvs2*rinvs2*rinvs2;

                const real lambdaFactor = (real(1.0)-lambda)*(real(1.0)-lambda);
                const real lambda_n     = pow(lambda,n);
                const real lambda_n_der = n*pow(lambda,n-1);

                const real denominator  = alpha*lambdaFactor+rinvs6;
                const real denominator2 = denominator*denominator;
                const real denominator3 = denominator2*denominator;

                const real ld = epsilon*(lambda_n_der/denominator2+(lambda_n*real(4.0)*alpha*(real(1.0)-lambda))/denominator3);

                return ld;
            }
        };
    }

}}}}
