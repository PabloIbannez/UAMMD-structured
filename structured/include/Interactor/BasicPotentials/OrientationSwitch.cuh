#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Utils/Maths/MatrixOperations.cuh"

#include "Interactor/BasicPotentials/StiffnessMask.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicPotentials{

    namespace OrientationSwitch{

        static inline __device__ real energy(const tensor3& A, const tensor3& B,
                                             const tensor3& R){

            tensor3 C = MatrixOperations::matmul(A.transpose(), B);
            const real frob = MatrixOperations::frobenius(C, R);

            return (real(1.0) + frob)/real(4.0);
        }

        static inline __device__ real3 torque(const tensor3& A, const tensor3& B,
                                              const tensor3& R){

            tensor3 dfrobdA = MatrixOperations::matmul(B, R.transpose());

            return  MatrixOperations::Pi(A, dfrobdA)/real(4.0);
        }

        static inline __device__ EnergyForceTorque energyForceTorque(const tensor3& A, const tensor3& B,
                                                                     const tensor3& R){
            EnergyForceTorque result;

            result.energy = energy(A, B, R);
            result.force  = make_real4(0.0);
            result.torque = make_real4(torque(A, B, R), 0.0);

            return result;
        }

        namespace Stiffness{

            static inline __device__ real energy(const tensor3& A, const tensor3& B,
                                                 const tensor3& R, const real& K){
                return StiffnessMask::stiffnessMask(OrientationSwitch::energy(A, B, R), K);
            }

            static inline __device__ real3 torque(const tensor3& A, const tensor3& B,
                                                  const tensor3& R, const real& K){

                const real e    = OrientationSwitch::energy(A, B, R);
                const real dUde = StiffnessMask::stiffnessMaskFirstDerivative(e, K);

                tensor3 dedA = MatrixOperations::matmul(B, R.transpose())/real(4.0);

                return dUde*MatrixOperations::Pi(A, dedA);
            }

            static inline __device__ EnergyForceTorque energyForceTorque(const tensor3& A, const tensor3& B,
                                                                         const tensor3& R, const real& K){

                const real4 force  = make_real4(0.0);

                const real e = OrientationSwitch::energy(A, B, R);

                const real dUde    = StiffnessMask::stiffnessMaskFirstDerivative(e, K);
                const tensor3 dedA = MatrixOperations::matmul(B, R.transpose())/real(4.0);

                EnergyForceTorque result;
                result.energy = StiffnessMask::stiffnessMask(e, K);
                result.force  = force;
                result.torque = make_real4(dUde*MatrixOperations::Pi(A, dedA), 0.0);

                return result;
            }
        }

    }


}}}}
