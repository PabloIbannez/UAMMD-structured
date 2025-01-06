#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Interactor/BasicPotentials/OrientationSwitch.cuh"
#include "Interactor/BasicPotentials/StiffnessMask.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

namespace RAP{

    // Being A, B rotation matrices, and R some target rotation matrix, the rotational alignment potential is defined as:
    // U = K*(real(3) - Frob(A^TB, R))/(real(4))

    static inline __device__ real energy(const tensor3& A, const tensor3& B,
                                         const tensor3& R, const real& K){

        const real swt = real(1.0)-OrientationSwitch::energy(A, B, R);

        return K*swt;
    }

    static inline __device__ real3 torque(const tensor3& A, const tensor3& B,
                                          const tensor3& R, const real& K){

        return -K*OrientationSwitch::torque(A, B, R);
    }

    static inline __device__ EnergyForceTorque energyForceTorque(const tensor3& A, const tensor3& B,
                                                                 const tensor3& R, const real& K){
        EnergyForceTorque result;

        result.energy = energy(A, B, R, K);
        result.force  = make_real4(0.0);
        result.torque = make_real4(torque(A, B, R, K), 0.0);

        return result;
    }

}}}}}
