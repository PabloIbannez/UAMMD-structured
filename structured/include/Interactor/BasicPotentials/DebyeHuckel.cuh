#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

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

}}}}
