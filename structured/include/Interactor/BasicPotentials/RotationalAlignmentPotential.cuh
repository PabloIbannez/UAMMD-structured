#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Utils/Maths/MatrixOperations.cuh"

#include "Interactor/BasicPotentials/StiffnessMask.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

namespace RAP{

    // Being A, B rotation matrices, and R some target rotation matrix, the rotational alignment potential is defined as:
    // U = (real(3) - Frob(A^TB, R))/(real(4))

    static inline __device__ real energy(const tensor3& A, const tensor3& B, const tensor3& R){
        tensor3 C = MatrixOperations::matmul(A.transpose(), B);
        const real frob = MatrixOperations::frobenius(C, R);
        return (real(3.0) - frob)/real(4.0);
    }

    static inline __device__ real3 torque(const tensor3& A, const tensor3& B, const tensor3& R){
        tensor3 dfrobdA    = MatrixOperations::matmul(B, R.transpose());
        return -MatrixOperations::Pi(A, dfrobdA)/real(4.0);
    }

    static inline __device__ EnergyForceTorque energyForceTorque(const tensor3& A, const tensor3& B, const tensor3& R){
        EnergyForceTorque result;
        result.energy = energy(A, B, R);
        result.force  = make_real4(0.0);
        result.torque = make_real4(torque(A, B, R), 0.0);
        return result;
    }

    namespace Stiffness{

        static inline __device__ real energy(const tensor3& A, const tensor3& B, const tensor3& R, const real& K){
            return StiffnessMask::stiffnessMask(RAP::energy(A, B, R), K);
        }

        static inline __device__ real3 torque(const tensor3& A, const tensor3& B, const tensor3& R, const real& K){
            const real e_rap    = RAP::energy(A, B, R);
            const real dUde_rap = StiffnessMask::stiffnessMaskFirstDerivative(e_rap, K);

            tensor3 de_rapdA = -MatrixOperations::matmul(B, R.transpose())/real(4.0);

            return dUde_rap*MatrixOperations::Pi(A, de_rapdA);
        }

        static inline __device__ EnergyForceTorque energyForceTorque(const tensor3& A, const tensor3& B, const tensor3& R, const real& K){
            const real4 force  = make_real4(0.0);

            const real e_rap = RAP::energy(A, B, R);

            const real dUde_rap    = StiffnessMask::stiffnessMaskFirstDerivative(e_rap, K);
            const tensor3 de_rapdA = -MatrixOperations::matmul(B, R.transpose())/real(4.0);

            EnergyForceTorque result;
            result.energy = StiffnessMask::stiffnessMask(e_rap, K);
            result.force  = force;
            result.torque = make_real4(dUde_rap*MatrixOperations::Pi(A, de_rapdA), 0.0);

            return result;
        }
    }

}}}}}
