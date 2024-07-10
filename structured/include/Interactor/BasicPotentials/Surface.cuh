#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

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

}}}}
