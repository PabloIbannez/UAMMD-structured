#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

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

}}}}
