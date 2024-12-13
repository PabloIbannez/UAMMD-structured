#pragma once

#include "utils/quaternion.cuh"
#include "Definitions/Computations.cuh"

#include "Utils/Maths/MatrixOperations.cuh"

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
        tensor3 denergydA = MatrixOperations::matmul(B, R.transpose());
        return -MatrixOperations::Pi(A, denergydA);
    }

}}}}}
