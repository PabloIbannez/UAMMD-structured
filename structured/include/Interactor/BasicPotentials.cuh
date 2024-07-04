#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    //Assumed force over i is computed

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
	    //TODO: Manage r = 0
            tensor3 H = tensor3(0.0);

	    const real invr2 = real(1.0)/r2;
            const real invr  = sqrt(invr2);

            H.xx =  K*(real(-1.0)+r0*invr*(real(1.0)-rij.x*rij.x*invr2));
            H.xy = -K*r0*invr*rij.x*rij.y*invr2;
            H.xz = -K*r0*invr*rij.x*rij.z*invr2;

            H.yx =  H.xy;
            H.yy =  K*(real(-1.0)+r0*invr*(real(1.0)-rij.y*rij.y*invr2));
            H.yz = -K*r0*invr*rij.y*rij.z*invr2;

            H.zx = H.xz;
            H.zy = H.yz;
            H.zz = K*(real(-1.0)+r0*invr*(real(1.0)-rij.z*rij.z*invr2));

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


    struct Morse{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){

            const real  r = sqrt(r2);
            const real dr = r-r0;
            const real factor = exp(-dr/D);

            const real fmod = real(2.0)*(e/D)*(real(1.0)-factor)*factor;

            return fmod*rij/r;
        }

        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){

            const real dr = sqrt(r2)-r0;
            const real factor = real(1.0)-exp(-dr/D);

            return e*(factor*factor-real(1.0));

        }
    };

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

    namespace WCA{

        struct Type1{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                //r2 > (sigma*2^(1/6))^2
                if(r2 > sigma*sigma*real(1.259921)) return make_real3(0);

                return LennardJones::Type1::force(rij,r2,epsilon,sigma);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 Stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    //r2 > (sigma*2^(1/6))^2
	    if(r2 > sigma*sigma*real(1.259921)) return real(0);

	    return LennardJones::Type1::hessian(rij,r2,epsilon,sigma);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                //r2 > (sigma*2^(1/6))^2
                if(r2 > sigma*sigma*real(1.259921)) return real(0);

                return LennardJones::Type1::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };

        struct Type2{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return make_real3(0);

                return LennardJones::Type2::force(rij,r2,epsilon,sigma);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }


	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    //r2 > (sigma*2^(1/6))^2
	    if(r2 > sigma*sigma*real(1.259921)) return real(0);

	    return LennardJones::Type2::hessian(rij,r2,epsilon,sigma);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return real(0);

                return LennardJones::Type2::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };

        struct Type3{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return make_real3(0);

                return LennardJones::Type3::force(rij,r2,epsilon,sigma);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    //r2 > (sigma*2^(1/6))^2
	    if(r2 > sigma*sigma*real(1.259921)) return real(0);

	    return LennardJones::Type3::hessian(rij,r2,epsilon,sigma);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return real(0);

                return LennardJones::Type3::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };
    }

    namespace GeneralLennardJones{

        struct Type1{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type1::force(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type1::force(rij,r2,fabs(epsilon),sigma);

            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    if(epsilon >= real(0.0)){
	      return WCA::Type1::hessian(rij,r2,epsilon,sigma);
	    }

	    return LennardJones::Type1::hessian(rij,r2,fabs(epsilon),sigma);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type1::energy(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type1::energy(rij,r2,fabs(epsilon),sigma);
            }
        };

        struct Type2{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type2::force(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type2::force(rij,r2,fabs(epsilon),sigma);

            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    if(epsilon >= real(0.0)){
	      return WCA::Type2::hessian(rij,r2,epsilon,sigma);
	    }

	    return LennardJones::Type2::hessian(rij,r2,fabs(epsilon),sigma);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type2::energy(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type2::energy(rij,r2,fabs(epsilon),sigma);
            }
        };

        struct Type3{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type3::force(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type3::force(rij,r2,fabs(epsilon),sigma);

            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& epsilon,const real& sigma){

	    if(epsilon >= real(0.0)){
	      return WCA::Type3::hessian(rij,r2,epsilon,sigma);
	    }

	    return LennardJones::Type3::hessian(rij,r2,fabs(epsilon),sigma);
	  }

	  //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type3::energy(rij,r2,epsilon,sigma);
                }

                return LennardJones::Type3::energy(rij,r2,fabs(epsilon),sigma);
            }
        };
    }

    namespace ModifiedLennardJones{

        struct Gaussian{

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){

                return WCA::Type2::force(rij,r2,epsilon,sigma)+
                       GaussianWell::force(rij,r2,epsilon,sigma,D);
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma,const real& D){
                return tensor3(0);
            }

	  //Energy
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						     const real& epsilon,const real& sigma,const real& D){
                return WCA::Type2::hessian(rij,r2,epsilon,sigma)+
                       GaussianWell::hessian(rij,r2,epsilon,sigma,D);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return WCA::Type2::energy(rij,r2,epsilon,sigma)+
                       GaussianWell::energy(rij,r2,epsilon,sigma,D);
            }

        };
    }

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

    struct NonPolar{

        //Force
        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){

            if(epsilon == real(0.0) ){
                return BasicPotentials::Steric::Steric::force<12>(rij,r2,zeroEnergy,sigma);
            }

            real3 f = BasicPotentials::LennardJones::Type1::force(rij,r2,fabs(epsilon),sigma);
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0) and r2>=r02){
                f=-f;
            }

            return f;
        }

        //Virial
        static inline __device__ real virial(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
            return real(0);
        }

        //Stress
        static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                const real& epsilon,const real& sigma,const real& zeroEnergy){
            return tensor3(0);
        }

      //Hessian
      static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
					       const real& epsilon,const real& sigma,const real& zeroEnergy){
	if(epsilon == real(0.0) ){
	  return BasicPotentials::Steric::Steric::hessian<12>(rij,r2,zeroEnergy,sigma);
	}

            tensor3 H = BasicPotentials::LennardJones::Type1::hessian(rij,r2,fabs(epsilon),sigma);
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0) and r2>=r02){
                H=-H;
            }

            return H;
      }

        //Energy
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){

            if(epsilon == real(0.0) ){
                return BasicPotentials::Steric::Steric::energy<12>(rij,r2,zeroEnergy,sigma);
            }

            real e  = BasicPotentials::LennardJones::Type1::energy(rij,r2,fabs(epsilon),sigma);

            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0)){

                if(r2<r02){
                    e+=real(2.0)*epsilon;
                } else {
                    e=-e;
                }
            }

            return e*real(0.5);
        }

    };

    namespace DebyeHuckel{

        struct DebyeHuckel{

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& dielectricConstant,const real& debyeLenght){

                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real efactor = ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("%f\n",efactor);

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real fmod = -efactor*exp(-r/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);

                //printf("%f\n",fmod);

                return fmod*rij;
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& dielectricConstant,const real& debyeLenght){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& ELECOEF,
                                                    const real& chgProduct,
                                                    const real& dielectricConstant,const real& debyeLenght){
                return tensor3(0);
            }

	  //hessian
            static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						     const real& ELECOEF,
						     const real& chgProduct,
						     const real& dielectricConstant,const real& debyeLenght){
	      const real      r  = sqrt(r2);
	      const real invr    = real(1.0)/r;
	      const real invr2   = real(1.0)/r2;

	      const real efactor = ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

	      real energyDerivative = -efactor*exp(-r/debyeLenght)*invr*(invr+real(1.0)/debyeLenght);

	      real energySecondDerivative = real(2.0)*efactor*exp(-r/debyeLenght)*invr*(invr2 +
											invr/debyeLenght+
											real(0.5)/(debyeLenght*debyeLenght));

	      return computeHessianRadialPotential(rij, invr, invr2, energyDerivative, energySecondDerivative);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& dielectricConstant,const real& debyeLenght){

	      const real      r  = sqrt(r2);

                const real efactor = ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real e = efactor*exp(-r/debyeLenght)/r;

                return e;
            }
        };

        struct DebyeHuckelSpheres{

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){

                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real efactor = (ELECOEF*chgProduct/dielectricConstant)/
                                     ((real(1.0)+radius1/debyeLenght)*(real(1.0)+radius2/debyeLenght)); //ELECOEF = 1/(4*pi*e_0)

                real fmod = -efactor*exp(-(r-radius1-radius2)/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);

                return fmod*rij;
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){
                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& ELECOEF,
                                                    const real& chgProduct,
                                                    const real& radius1, const real& radius2,
                                                    const real& dielectricConstant,const real& debyeLenght){
                return tensor3(0);
            }

	  //Hessian
	  static inline __device__ tensor3 hessian(const real3& rij, const real& r2,
						   const real& ELECOEF,
						   const real& chgProduct,
						   const real& radius1, const real& radius2,
						   const real& dielectricConstant,const real& debyeLenght){

	    const real      r  = sqrt(r2);
	    const real invr    = real(1.0)/r;
	    const real invr2   = real(1.0)/r2;

	    const real dlenght2 = debyeLenght*debyeLenght;

	    const real efactor = (ELECOEF*chgProduct/dielectricConstant)/
	                         ((real(1.0)+radius1/debyeLenght)*(real(1.0)+radius2/debyeLenght)); //ELECOEF = 1/(4*pi*e_0)


	    real energyDerivative = -efactor*(debyeLenght + r)*exp(-(r-radius1-radius2)/debyeLenght)*invr2/debyeLenght;

	    real energySecondDerivative = efactor*(real(2.0)*dlenght2 + real(2.0)*debyeLenght*r + r2)*exp(-(r-radius1-radius2)/debyeLenght)*invr*invr2/dlenght2;


	    return computeHessianRadialPotential(rij, invr, invr2, energyDerivative, energySecondDerivative);
	  }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){

                const real      r  = sqrt(r2);

                const real efactor = (ELECOEF*chgProduct/dielectricConstant)/
                                     ((real(1.0)+radius1/debyeLenght)*(real(1.0)+radius2/debyeLenght)); //ELECOEF = 1/(4*pi*e_0)

                real e = efactor*exp(-(r-radius1-radius2)/debyeLenght)/r;

                return e;
            }
        };

        inline __device__ real distanceDependentDielectric(const real& r,
                                                           const real& A, const real& B,
                                                           const real& k, const real lmd){
            return A+B/(real(1.0)+k*exp(-lmd*B*r));
        }

        inline __device__ real distanceDependentDielectricDerivate(const real& r,
                                                                   const real& A, const real& B,
                                                                   const real& k, const real lmd){

            real den = exp(lmd*B*r)+k;
            return B*B*k*lmd*exp(lmd*B*r)/(den*den);
        }

        struct DebyeHuckelDistanceDependentDielectric{

            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){

                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real dielectricConstant = distanceDependentDielectric(r,A,B,k,lmd);

                const real efactor = ELECOEF*chgProduct/dielectricConstant;

                real fmod = -efactor*exp(-r/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);
                     fmod = fmod - (efactor*exp(-r/debyeLenght)/r)*distanceDependentDielectricDerivate(r,A,B,k,lmd)/(r*dielectricConstant);

                return fmod*rij;
            }

            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){

                return real(0);
            }

            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& ELECOEF,
                                                    const real& chgProduct,const real& debyeLenght,
                                                    const real& A, const real& B,
                                                    const real& k, const real lmd){

                return tensor3(0);
            }

            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& ELECOEF,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){

                const real      r  = sqrt(r2);

                const real dielectricConstant = distanceDependentDielectric(r,A,B,k,lmd);

                const real efactor = ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real e = efactor*exp(-r/debyeLenght)/r;

                return e;
            }
        };

    }

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

    namespace Surface{

        struct Parabola{

            //Force
            static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                 const real& epsilon){

                real dz=real(0.0);
                if(pos.z <= surfacePos){
                    dz = surfacePos - pos.z;
                }

                return make_real3(0.0,0.0,epsilon*dz);
            }

            //Energy
            static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                 const real& epsilon){

                real dz2=real(0.0);
                if(pos.z <= surfacePos){
                    dz2 = surfacePos - pos.z;
                    dz2 = dz2*dz2;
                }

                return real(0.5)*epsilon*dz2;
            }
        };

        struct HarmonicWell{

            //Force
            static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                 const real& epsilon,const real& sigma){

                const real z = pos.z;
                if( (z < (-sigma + surfacePos)) or (z > (sigma + surfacePos))){
                    return make_real3(0.0,0.0,0.0);
                }

                return make_real3(0.0,0.0,epsilon*(surfacePos-z)/sigma/sigma);
            }

            //Energy
            static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                 const real& epsilon,const real& sigma){

                const real z = pos.z;
                if( (z < (-sigma + surfacePos)) or (z > (sigma + surfacePos))){
                    return real(0.0);
                }

                const real d = (surfacePos-z)/sigma;
                const real e = real(0.5)*epsilon*(d*d-real(1.0));

                return e;
            }
        };

        struct GaussianWell{

            static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                 const real& e,const real& D){

                real dz = pos.z - surfacePos;

                return -e*exp(-dz*dz/(real(2.0)*D));

            }

            static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                 const real& e,const real& D){

                real dz = pos.z - surfacePos;
                real Fz = -dz * e * exp(-dz * dz / (real(2.0) * D)) / D;

                return make_real3(real(0.0), real(0.0), Fz);
            }


        };

        namespace LennardJones{

            struct Type1{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    const real invdz2   = real(1.0)/(dz*dz);
                    const real sinvdz2  = sigma*sigma*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real fmod = real(24.0)*epsilon*(real(2.0)*sinvdz12-sinvdz6)*invdz2;

                    return make_real3(0.0,0.0,fmod*(pos.z-surfacePos));
                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    const real invdz2   = real(1.0)/(dz*dz);
                    const real sinvdz2  = sigma*sigma*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    return real(4.0)*epsilon*(sinvdz12-sinvdz6);
                }
            };

            struct Type2{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    const real invdz2   = real(1.0)/(dz*dz);
                    const real sinvdz2  = sigma*sigma*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real fmod = real(12.0)*epsilon*(sinvdz12-sinvdz6)*invdz2;

                    return make_real3(0.0,0.0,fmod*(pos.z-surfacePos));
                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    const real invdz2   = real(1.0)/(dz*dz);
                    const real sinvdz2  = sigma*sigma*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    return epsilon*(sinvdz12-real(2.0)*sinvdz6);
                }
            };
        }

        namespace WCA{

            struct Type1{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    if(dz > sigma*real(1.122462)) return make_real3(0);

                    return  LennardJones::Type1::force(pos,surfacePos,epsilon,sigma);

                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    if(dz > sigma*real(1.122462)) return real(0);

                    return  LennardJones::Type1::energy(pos,surfacePos,epsilon,sigma) + epsilon;

                }
            };

            struct Type2{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    if(dz > sigma) return make_real3(0);

                    return  LennardJones::Type2::force(pos,surfacePos,epsilon,sigma);

                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    const real dz = fabs(pos.z - surfacePos);

                    if(dz > sigma) return real(0);

                    return  LennardJones::Type2::energy(pos,surfacePos,epsilon,sigma) + epsilon;

                }
            };



        }

        namespace GeneralLennardJones{

            struct Type1{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    if(epsilon >= real(0.0)){
                        return WCA::Type1::force(pos,surfacePos,epsilon,sigma);
                    }

                    return LennardJones::Type1::force(pos,surfacePos,fabs(epsilon),sigma);

                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    if(epsilon >= real(0.0)){
                        return WCA::Type1::energy(pos,surfacePos,epsilon,sigma);
                    }

                    return LennardJones::Type1::energy(pos,surfacePos,fabs(epsilon),sigma);

                }
            };

            struct Type2{

                //Force
                static inline __device__ real3 force(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    if(epsilon >= real(0.0)){
                        return WCA::Type2::force(pos,surfacePos,epsilon,sigma);
                    }

                    return LennardJones::Type2::force(pos,surfacePos,fabs(epsilon),sigma);

                }

                //Energy
                static inline __device__ real energy(const real3& pos, const real& surfacePos,
                                                     const real& epsilon,const real& sigma){

                    if(epsilon >= real(0.0)){
                        return WCA::Type2::energy(pos,surfacePos,epsilon,sigma);
                    }

                    return LennardJones::Type2::energy(pos,surfacePos,fabs(epsilon),sigma);

                }
            };
        }
    }

    namespace Tip{

        struct GenericSphericalTip{

            //Force
            static inline __device__ real3 force(const real3& pos, const real3& tipPos,
                                                 const real& A ,const real& B,
                                                 const real& tipRadius,const real& epsilon,const real& sigma,
                                                 Box box){

                const real3 dr = box.apply_pbc(tipPos - pos);

                const real r  = sqrt(dot(dr,dr));
                      real r2 = r-tipRadius;
                           r2 = r2*r2;

                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -epsilon*real(6.0)*(real(2.0)*A*sinvr12+B*sinvr6)/abs(r-tipRadius);

                real3 f = fmod*dr/r;

                return f;
            }

            //Energy
            static inline __device__ real energy(const real3& pos, const real3& tipPos,
                                                 const real& A ,const real& B,
                                                 const real& tipRadius,const real& epsilon,const real& sigma,
                                                 Box box){

                const real3 dr = box.apply_pbc(tipPos - pos);

                const real r  = sqrt(dot(dr,dr));
                      real r2 = r-tipRadius;
                           r2 = r2*r2;

                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(A*sinvr12+B*sinvr6);

                return e;

            }
        };
    }

    namespace Helix{

        static inline __device__ real3 ex_optimal(const Quat& qs,
                                                  const tensor3 R_H){
            const real3 ex_opt = rotateVector(qs,make_real3(R_H.xx,R_H.yx,R_H.zx));
            return ex_opt;
        }

        static inline __device__ real3 ey_optimal(const Quat& qs,
                                                  const tensor3 R_H){
            const real3 ey_opt = rotateVector(qs,make_real3(R_H.xy,R_H.yy,R_H.zy));
            return ey_opt;
        }

        static inline __device__ real cos_theta(const Quat& qs, const Quat& qe,
                                                const tensor3 R_H){

            const real3 ex     = qe.getVx();
            const real3 ex_opt = ex_optimal(qs,R_H);

            real ctheta = dot(ex,ex_opt);

            ctheta = fminf(ctheta, real(1.0));
            ctheta = fmaxf(ctheta,-real(1.0));

            return ctheta;
        }

        static inline __device__ real cos_phi(const Quat& qs, const Quat& qe,
                                              const tensor3 R_H){

            const real3 ey     = qe.getVy();
            const real3 ey_opt = ey_optimal(qs,R_H);

            real cphi = dot(ey,ey_opt);

            cphi = fminf(cphi, real(1.0));
            cphi = fmaxf(cphi,-real(1.0));

            return cphi;
        }

        struct exponential{

            struct potential{

                struct cosParameters{
                    real K;
                };

                static inline __device__ real cosEnergy(const Quat& qs, const Quat& qe,
                                                        const real& cosine,const cosParameters& params){

                    const real K = params.K;
                    if(K == real(0.0)) return real(0.0);

                    const real A   = expf(-real(2.0)*K);
                    const real den = real(1.0)-A;        // Denominator

                    const real e = (expf(-K*(real(1.0)-cosine))-A)/den;

                    return e;
                }

                static inline __device__ real cosEnergyDerivative(const Quat& qs, const Quat& qe,
                                                                  const real& cosine,const cosParameters& params){

                    const real K = params.K;
                    if(K == real(0.0)) return real(0.0);

                    const real A   = expf(-real(2.0)*K);
                    const real den = real(1.0)-A;        // Denominator
                    const real dU_dcos_theta = K*expf(-K*(real(1.0)-cosine))/den;

                    return dU_dcos_theta;
                }

                struct dstParameters{
                    real rc;
                    real K;
                };

                static inline __device__ real dstEnergy(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchExponential::energy(dr,r2,real(-1.0),params.rc,params.K);
                }

                static inline __device__ real3 dstForce(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchExponential::force(dr,r2,real(-1.0),params.rc,params.K);
                }

            };

            struct params{
                typename potential::dstParameters dstParams;
                typename potential::cosParameters thetaParams;
                typename potential::cosParameters phiParams;
            };

            static params readParams(DataEntry& data, bool loadDst = false){
                params p;

                p.thetaParams.K = data.getParameter<real>("Ka");
                p.phiParams.K   = data.getParameter<real>("Kd");

                System::log<System::MESSAGE>("[Helix] Ka = %f", p.thetaParams.K);
                System::log<System::MESSAGE>("[Helix] Kd = %f", p.phiParams.K);

                if(loadDst){
                    p.dstParams.rc = data.getParameter<real>("rc");
                    p.dstParams.K  = data.getParameter<real>("Kb");

                    System::log<System::MESSAGE>("[Helix] rc = %f", p.dstParams.rc);
                    System::log<System::MESSAGE>("[Helix] Kb = %f", p.dstParams.K);
                }

                return p;
            }

            template <typename T>
            static params readParamsMap(std::map<std::string,T>& info,std::string suffix = ""){
                params p;

                p.dstParams.rc = real(info.at("rc"+suffix));
                p.dstParams.K  = real(info.at("Kb"+suffix));

                p.thetaParams.K = real(info.at("Ka"+suffix));
                p.phiParams.K   = real(info.at("Kd"+suffix));

                System::log<System::MESSAGE>("[Helix] rc%s = %f", suffix.c_str(), p.dstParams.rc);
                System::log<System::MESSAGE>("[Helix] Kb%s = %f", suffix.c_str(), p.dstParams.K);

                System::log<System::MESSAGE>("[Helix] Ka%s = %f", suffix.c_str(), p.thetaParams.K);
                System::log<System::MESSAGE>("[Helix] Kd%s = %f", suffix.c_str(), p.phiParams.K);

                return p;
            }
        };

        struct cosine{

            struct potential{

                struct cosParameters{
                    real cos_angle_start;
                    real cos_angle_end;
                };

                static inline __device__ real cosEnergy(const Quat& qs, const Quat& qe,
                                                        const real& cosine,const cosParameters& params){

                    const real cos_angle_start = params.cos_angle_start;
                    const real cos_angle_end = params.cos_angle_end;

                    real swt;
                    if(cosine > cos_angle_start){
                        swt = real(-1.0);
                    } else if(cosine < cos_angle_end){
                        swt = real(1.0);
                    } else {
                        real norm = (cosine-cos_angle_start)/(cos_angle_end-cos_angle_start);

                        swt = cosf(norm*real(M_PI));

                        swt = fminf(swt, real(1.0));
                        swt = fmaxf(swt,-real(1.0));

                        swt = -swt;
                    }

                    swt = (real(1.0)-swt)*real(0.5);

                    return swt;
                }

                static inline __device__ real cosEnergyDerivative(const Quat& qs, const Quat& qe,
                                                                  const real& cosine,const cosParameters& params){

                    const real cos_angle_start = params.cos_angle_start;
                    const real cos_angle_end = params.cos_angle_end;

                    real dswt_dcos;

                    if(cosine > cos_angle_start){
                        dswt_dcos = real(0.0);
                    } else if(cosine < cos_angle_end){
                        dswt_dcos = real(0.0);
                    } else {
                        real norm = (cosine-cos_angle_start)/(cos_angle_end-cos_angle_start);

                        dswt_dcos = sinf(norm*real(M_PI));

                        dswt_dcos = fminf(dswt_dcos, real(1.0));
                        dswt_dcos = fmaxf(dswt_dcos,-real(1.0));

                        dswt_dcos = dswt_dcos*real(M_PI)/(cos_angle_end-cos_angle_start);
                    }

                    dswt_dcos = -dswt_dcos*real(0.5);

                    return dswt_dcos;
                }

                struct dstParameters{
                    real r_start;
                    real rc;
                };

                static inline __device__ real dstEnergy(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchCosine::energy(dr,r2,real(-1.0),params.r_start,params.rc);
                }

                static inline __device__ real3 dstForce(const real3 dr,
                                                        const real  r2,
                                                        const dstParameters& params){
                    return DistanceSwitchCosine::force(dr,r2,real(-1.0),params.r_start,params.rc);
                }

            };

            struct params{
                typename potential::dstParameters dstParams;
                typename potential::cosParameters thetaParams;
                typename potential::cosParameters phiParams;
            };

            static params readParams(DataEntry& data, bool loadDst = false){
                params p;

                p.thetaParams.cos_angle_start = data.getParameter<real>("theta_start");
                p.thetaParams.cos_angle_end   = data.getParameter<real>("theta_end");

                p.phiParams.cos_angle_start = data.getParameter<real>("phi_start");
                p.phiParams.cos_angle_end   = data.getParameter<real>("phi_end");

                System::log<System::MESSAGE>("[Helix] theta_start = %f", p.thetaParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] theta_end   = %f", p.thetaParams.cos_angle_end);

                System::log<System::MESSAGE>("[Helix] phi_start = %f", p.phiParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] phi_end   = %f", p.phiParams.cos_angle_end);

                p.thetaParams.cos_angle_start = cos(p.thetaParams.cos_angle_start);
                p.thetaParams.cos_angle_start = min(p.thetaParams.cos_angle_start, real(1.0));
                p.thetaParams.cos_angle_start = max(p.thetaParams.cos_angle_start,-real(1.0));

                p.thetaParams.cos_angle_end   = cos(p.thetaParams.cos_angle_end);
                p.thetaParams.cos_angle_end   = min(p.thetaParams.cos_angle_end, real(1.0));
                p.thetaParams.cos_angle_end   = max(p.thetaParams.cos_angle_end,-real(1.0));

                p.phiParams.cos_angle_start = cos(p.phiParams.cos_angle_start);
                p.phiParams.cos_angle_start = min(p.phiParams.cos_angle_start, real(1.0));
                p.phiParams.cos_angle_start = max(p.phiParams.cos_angle_start,-real(1.0));

                p.phiParams.cos_angle_end   = cos(p.phiParams.cos_angle_end);
                p.phiParams.cos_angle_end   = min(p.phiParams.cos_angle_end, real(1.0));
                p.phiParams.cos_angle_end   = max(p.phiParams.cos_angle_end,-real(1.0));

                if(loadDst){
                    p.dstParams.r_start = data.getParameter<real>("r_start");
                    p.dstParams.rc      = data.getParameter<real>("rc");

                    System::log<System::MESSAGE>("[Helix] r_start = %f", p.dstParams.r_start);
                    System::log<System::MESSAGE>("[Helix] rc      = %f", p.dstParams.rc);
                }

                return p;
            }

            template <typename T>
            static params readParamsMap(std::map<std::string,T>& info,std::string suffix = ""){
                params p;

                p.dstParams.rc = real(info.at("rc"+suffix));
                p.dstParams.r_start = real(info.at("r_start"+suffix));

                p.thetaParams.cos_angle_start = real(info.at("theta_start"+suffix));
                p.thetaParams.cos_angle_end   = real(info.at("theta_end"+suffix));

                p.phiParams.cos_angle_start = real(info.at("phi_start"+suffix));
                p.phiParams.cos_angle_end   = real(info.at("phi_end"+suffix));

                System::log<System::MESSAGE>("[Helix] rc%s = %f", suffix.c_str(), p.dstParams.rc);
                System::log<System::MESSAGE>("[Helix] r_start%s = %f", suffix.c_str(), p.dstParams.r_start);

                System::log<System::MESSAGE>("[Helix] theta_start%s = %f", suffix.c_str(), p.thetaParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] theta_end%s   = %f", suffix.c_str(), p.thetaParams.cos_angle_end);

                System::log<System::MESSAGE>("[Helix] phi_start%s = %f", suffix.c_str(), p.phiParams.cos_angle_start);
                System::log<System::MESSAGE>("[Helix] phi_end%s   = %f", suffix.c_str(), p.phiParams.cos_angle_end);

                p.thetaParams.cos_angle_start = cos(p.thetaParams.cos_angle_start);
                p.thetaParams.cos_angle_start = min(p.thetaParams.cos_angle_start, real(1.0));
                p.thetaParams.cos_angle_start = max(p.thetaParams.cos_angle_start,-real(1.0));

                p.thetaParams.cos_angle_end   = cos(p.thetaParams.cos_angle_end);
                p.thetaParams.cos_angle_end   = min(p.thetaParams.cos_angle_end, real(1.0));
                p.thetaParams.cos_angle_end   = max(p.thetaParams.cos_angle_end,-real(1.0));

                p.phiParams.cos_angle_start = cos(p.phiParams.cos_angle_start);
                p.phiParams.cos_angle_start = min(p.phiParams.cos_angle_start, real(1.0));
                p.phiParams.cos_angle_start = max(p.phiParams.cos_angle_start,-real(1.0));

                p.phiParams.cos_angle_end   = cos(p.phiParams.cos_angle_end);
                p.phiParams.cos_angle_end   = min(p.phiParams.cos_angle_end, real(1.0));
                p.phiParams.cos_angle_end   = max(p.phiParams.cos_angle_end,-real(1.0));

                return p;
            }
        };

        template <typename potential>
        static inline __device__ real energyTheta(const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H, const typename potential::cosParameters& params){



            real ctheta = cos_theta(qs,qe,R_H);
            return potential::cosEnergy(qs,qe,ctheta,params);

        }

        template <typename potential>
        static inline __device__ real energyThetaForwardBackward(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                                 const typename potential::cosParameters& params){

            const real energyForward  = energyTheta<potential>(qs,qe,R_H,params);
            const real energyBackward = energyTheta<potential>(qe,qs,R_H.transpose(),params);

            return (energyForward+energyBackward)*real(0.5);
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativeTheta(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                          const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dsx
            // Second row: dU/dsy
            // Third row:  dU/dsz

            const real ctheta = cos_theta(qs,qe,R_H);
            const real dU_dcos_theta = potential::cosEnergyDerivative(qs,qe,ctheta,params);

            const real3 ex = qe.getVx();
            const real3 dcos_theta_dsx = ex*R_H.xx;
            const real3 dcos_theta_dsy = ex*R_H.yx;
            const real3 dcos_theta_dsz = ex*R_H.zx;

            const real3 dU_dsx = dU_dcos_theta*dcos_theta_dsx;
            const real3 dU_dsy = dU_dcos_theta*dcos_theta_dsy;
            const real3 dU_dsz = dU_dcos_theta*dcos_theta_dsz;

            tensor3 dU_dS;

            dU_dS.xx = dU_dsx.x;
            dU_dS.xy = dU_dsx.y;
            dU_dS.xz = dU_dsx.z;

            dU_dS.yx = dU_dsy.x;
            dU_dS.yy = dU_dsy.y;
            dU_dS.yz = dU_dsy.z;

            dU_dS.zx = dU_dsz.x;
            dU_dS.zy = dU_dsz.y;
            dU_dS.zz = dU_dsz.z;

            return dU_dS;
        }

        template<typename potential>
        static inline __device__ tensor3 eDerivativeTheta(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                          const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dex
            // Second row: dU/dey
            // Third row:  dU/dez

            const real ctheta = cos_theta(qs,qe,R_H);
            const real dU_dcos_theta = potential::cosEnergyDerivative(qs,qe,ctheta,params);

            const real3 dcos_theta_dex = ex_optimal(qs,R_H);
            //const real3 dcos_theta_dey = make_real3(0.0);
            //const real3 dcos_theta_dez = make_real3(0.0);

            const real3 dU_dex = dU_dcos_theta*dcos_theta_dex;

            tensor3 dU_dE;

            dU_dE.xx = dU_dex.x;
            dU_dE.xy = dU_dex.y;
            dU_dE.xz = dU_dex.z;

            dU_dE.yx = real(0.0);
            dU_dE.yy = real(0.0);
            dU_dE.yz = real(0.0);

            dU_dE.zx = real(0.0);
            dU_dE.zy = real(0.0);
            dU_dE.zz = real(0.0);

            return dU_dE;
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativeThetaForwardBackward(const Quat& qs, const Quat& qe,
                                                                         const tensor3 R_H,
                                                                         const typename potential::cosParameters& params){

            const tensor3 dU_ds_forward  = sDerivativeTheta<potential>(qs,qe,R_H,params);
            const tensor3 dU_ds_backward = eDerivativeTheta<potential>(qe,qs,R_H.transpose(),params);

            const tensor3 dU_ds = (dU_ds_forward + dU_ds_backward)*real(0.5);

            return dU_ds;
        }

        // PHI

        template <typename potential>
        static inline __device__ real energyPhi(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                const typename potential::cosParameters& params){

            real cphi = cos_phi(qs,qe,R_H);
            return potential::cosEnergy(qs,qe,cphi,params);
        }

        template <typename potential>
        static inline __device__ real energyPhiForwardBackward(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                               const typename potential::cosParameters& params){

            const real energyForward  = energyPhi<potential>(qs,qe,R_H,params);
            const real energyBackward = energyPhi<potential>(qe,qs,R_H.transpose(),params);

            return (energyForward+energyBackward)*real(0.5);
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativePhi(const Quat& qs, const Quat& qe,
                                                        const tensor3 R_H,
                                                        const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dsx
            // Second row: dU/dsy
            // Third row:  dU/dsz

            const real cphi = cos_phi(qs,qe,R_H);
            const real dU_dcos_phi = potential::cosEnergyDerivative(qs,qe,cphi,params);

            const real3 ey = qe.getVy();
            const real3 dcos_phi_dsx = ey*R_H.xy;
            const real3 dcos_phi_dsy = ey*R_H.yy;
            const real3 dcos_phi_dsz = ey*R_H.zy;

            const real3 dU_dsx = dU_dcos_phi*dcos_phi_dsx;
            const real3 dU_dsy = dU_dcos_phi*dcos_phi_dsy;
            const real3 dU_dsz = dU_dcos_phi*dcos_phi_dsz;

            tensor3 dU_dS;

            dU_dS.xx = dU_dsx.x;
            dU_dS.xy = dU_dsx.y;
            dU_dS.xz = dU_dsx.z;

            dU_dS.yx = dU_dsy.x;
            dU_dS.yy = dU_dsy.y;
            dU_dS.yz = dU_dsy.z;

            dU_dS.zx = dU_dsz.x;
            dU_dS.zy = dU_dsz.y;
            dU_dS.zz = dU_dsz.z;

            return dU_dS;
        }

        template <typename potential>
        static inline __device__ tensor3 eDerivativePhi(const Quat& qs, const Quat& qe, const tensor3 R_H,
                                                        const typename potential::cosParameters& params){

            // It returns tensor3:
            // First row:  dU/dex
            // Second row: dU/dey
            // Third row:  dU/dez

            const real cphi = cos_phi(qs,qe,R_H);
            const real dU_dcos_phi = potential::cosEnergyDerivative(qs,qe,cphi,params);

            //const real3 dcos_phi_dex = make_real3(0.0);
            const real3 dcos_phi_dey   = ey_optimal(qs,R_H);
            //const real3 dcos_phi_dez = make_real3(0.0);

            const real3 dU_dey = dU_dcos_phi*dcos_phi_dey;

            tensor3 dU_dE;

            dU_dE.xx = real(0.0);
            dU_dE.xy = real(0.0);
            dU_dE.xz = real(0.0);

            dU_dE.yx = dU_dey.x;
            dU_dE.yy = dU_dey.y;
            dU_dE.yz = dU_dey.z;

            dU_dE.zx = real(0.0);
            dU_dE.zy = real(0.0);
            dU_dE.zz = real(0.0);

            return dU_dE;
        }

        template <typename potential>
        static inline __device__ tensor3 sDerivativePhiForwardBackward(const Quat& qs, const Quat& qe,
                                                                       const tensor3 R_H,
                                                                       const typename potential::cosParameters& params){


            const tensor3 dU_ds_forward  = sDerivativePhi<potential>(qs,qe,R_H,params);
            const tensor3 dU_ds_backward = eDerivativePhi<potential>(qe,qs,R_H.transpose(),params);

            const tensor3 dU_ds = (dU_ds_forward + dU_ds_backward)*real(0.5);

            return dU_ds;
        }

        // Combined theta and phi

        template <typename potential>
        static inline __device__ real orientationEnergy(const Quat& qs, const Quat& qe,
                                                        const tensor3 R_H,
                                                        const typename potential::cosParameters& paramsTheta,
                                                        const typename potential::cosParameters& paramsPhi){
            //The orientation energy is given by:
            // U = UthetaFB*UphiFB

            return energyThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*energyPhiForwardBackward<potential>(qs,qe,R_H,paramsPhi);
        }

        //For orientation only no forces are present

        template <typename potential>
        static inline __device__ real3 orientationTorque(const Quat& qs, const Quat& qe,
                                                         const tensor3 R_H,
                                                         const typename potential::cosParameters& paramsTheta,
                                                         const typename potential::cosParameters& paramsPhi){

            //The orientation energy is given by:
            // U = UthetaFB*UphiFB
            // This function returns the total torque over the particle s


            // First we compute derivates of the potential energy respect s
            // dU/ds = dUthetaFB/ds*UphiFB + UthetaFB*dUphiFB/ds

            const tensor3 dU_ds = sDerivativeThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*energyPhiForwardBackward<potential>(qs,qe,R_H,paramsPhi) +
                                  energyThetaForwardBackward<potential>(qs,qe,R_H,paramsTheta)*sDerivativePhiForwardBackward<potential>(qs,qe,R_H,paramsPhi);

            // Then we compute the torque as:
            // T = -sx X dU/dsx - sy X dU/dsy - sz X dU/dsz

            real3 T = make_real3(0.0);

            const real3 dUdsx = make_real3(dU_ds.xx,dU_ds.xy,dU_ds.xz);
            const real3 dUdsy = make_real3(dU_ds.yx,dU_ds.yy,dU_ds.yz);
            const real3 dUdsz = make_real3(dU_ds.zx,dU_ds.zy,dU_ds.zz);

            const real3 sx = qs.getVx();
            const real3 sy = qs.getVy();
            const real3 sz = qs.getVz();

            T += -cross(sx,dUdsx);
            T += -cross(sy,dUdsy);
            T += -cross(sz,dUdsz);

            return T;

        }

        namespace Fixed {

            template <typename potential>
            static inline __device__ real energy(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const real& Kb,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){

                // dr is expected to be dr = e - s

                //The fixed energy is given by:
                // U = -E*UthetaFB*UphiFB+0.5*Kb*(r)^2
                // Where r is the distance between the points s and e
                // r^2 = dot(dr,dr)

                const real r2 = dot(dr,dr);

                return -E*orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi)
                        + Harmonic::energy(dr,r2,Kb,real(0.0));
            }

            static inline __device__ real3 force(const real3& dr,
                                                 const real& Kb){

                // Only the distance dependent of the fixed energy contributes to the force

                const real r2 = dot(dr,dr);

                if (r2 < real(1e-6)) {
                    return make_real3(0.0);
                }

                return Harmonic::force(dr,r2,Kb,real(0.0));
            }

            template <typename potential>
            static inline __device__ real3 torque(const real3& dr,
                                                  const real3& lps, // Local positions of s
                                                  const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H,
                                                  const real& Kb,
                                                  const typename potential::potential::cosParameters& paramsTheta,
                                                  const typename potential::potential::cosParameters& paramsPhi,
                                                  const real& E){

                // lps is a vector from the particle center to the point s

                // We first compute the torque due to the orientation

                const real3 T  = -E*orientationTorque<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
                const real3 fs =    force(dr,Kb);

                return T + cross(lps,fs);
            }

        } // namespace Fixed

        namespace Dynamic {

            template <typename potential>
            static inline __device__ real energy(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const typename potential::potential::dstParameters& paramsDst,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){
                // dr is expected to be dr = e - s

                //The dynamic energy is given by:
                // U = -E*Ubond*UthetaFB*UphiFB
                const real r2 = dot(dr,dr);

                const real Ub = -E*potential::potential::dstEnergy(dr,r2,paramsDst);
                const real Uo = orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);

                return Ub*Uo;
            }

            template <typename potential>
            static inline __device__ real3 force(const real3& dr,
                                                 const Quat& qs, const Quat& qe,
                                                 const tensor3 R_H,
                                                 const typename potential::potential::dstParameters& paramsDst,
                                                 const typename potential::potential::cosParameters& paramsTheta,
                                                 const typename potential::potential::cosParameters& paramsPhi,
                                                 const real& E){

                // dr is expected to be dr = e - s

                //The dynamic energy is given by:
                // U = -E*Ubond*UthetaFB*UphiFB
                // Only the distance dependent (Ubond) of the dynamic energy contributes to the force
                // fs = -dU/dr = -E*(-dUbond/dr)*UthetaFB*UphiFB = -E*dstSwitch::force*orientationEnergy

                const real r2 = dot(dr,dr);

                return -E*potential::potential::dstForce(dr,r2,paramsDst)*orientationEnergy<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
            }

            template <typename potential>
            static inline __device__ real3 torque(const real3& dr,
                                                  const real3& lps, // Local positions of s
                                                  const Quat& qs, const Quat& qe,
                                                  const tensor3 R_H,
                                                  const typename potential::potential::dstParameters& paramsDst,
                                                  const typename potential::potential::cosParameters& paramsTheta,
                                                  const typename potential::potential::cosParameters& paramsPhi,
                                                  const real& E){

                // lps is a vector from the particle center to the point s

                // U = -E*Ubond*UthetaFB*UphiFB

                // dU/ds = -EUbond*(dUthetaFB/ds*UphiFB + UthetaFB*dUphiFB/ds)

                const real r2 = dot(dr,dr);

                const real Ubond = -E*potential::potential::dstEnergy(dr,r2,paramsDst);

                const real3 T  =  Ubond*orientationTorque<potential::potential>(qs,qe,R_H,paramsTheta,paramsPhi);
                const real3 fs =  force<potential>(dr,qs,qe,R_H,paramsDst,paramsTheta,paramsPhi,E);

                return T + cross(lps,fs);
            }
        }
    }


}}}}

