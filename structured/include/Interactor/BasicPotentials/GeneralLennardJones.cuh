#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/LennardJones.cuh"
#include "Interactor/BasicPotentials/WCA.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

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

}}}}
