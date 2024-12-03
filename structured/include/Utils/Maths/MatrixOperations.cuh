#pragma once

#include "uammd.cuh"

#include "utils/quaternion.cuh"
#include "utils/tensor.cuh"

namespace uammd{
namespace structured{
namespace MatrixOperations{

    // Quaternion to matrix

    inline __host__ __device__ tensor3 quat2mat(const Quat& q){
        const real3 vx = q.getVx();
        const real3 vy = q.getVy();
        const real3 vz = q.getVz();

        tensor3 mat;

        mat.xx = vx.x;
        mat.yx = vx.y;
        mat.zx = vx.z;

        mat.xy = vy.x;
        mat.yy = vy.y;
        mat.zy = vy.z;

        mat.xz = vz.x;
        mat.yz = vz.y;
        mat.zz = vz.z;

        return mat;
    }

    inline __host__ __device__ tensor3 quat2mat(const real4& q){
        return quat2mat(Quat(q));
    }

    // Matrix multiplication

    inline __host__ __device__ tensor3 matmul(const tensor3& A, const tensor3& B){

        tensor3 R;

        R.xx = A.xx*B.xx + A.xy*B.yx + A.xz*B.zx;
        R.xy = A.xx*B.xy + A.xy*B.yy + A.xz*B.zy;
        R.xz = A.xx*B.xz + A.xy*B.yz + A.xz*B.zz;

        R.yx = A.yx*B.xx + A.yy*B.yx + A.yz*B.zx;
        R.yy = A.yx*B.xy + A.yy*B.yy + A.yz*B.zy;
        R.yz = A.yx*B.xz + A.yy*B.yz + A.yz*B.zz;

        R.zx = A.zx*B.xx + A.zy*B.yx + A.zz*B.zx;
        R.zy = A.zx*B.xy + A.zy*B.yy + A.zz*B.zy;
        R.zz = A.zx*B.xz + A.zy*B.yz + A.zz*B.zz;

        return R;
    }

    // Frobenius product

    inline __host__ __device__ real frobenius(const tensor3& A, const tensor3& B){
        // Frob(A,B) = tr(A^T B)

        const real xx = A.xx*B.xx + A.yx*B.yx + A.zx*B.zx;
        const real yy = A.xy*B.xy + A.yy*B.yy + A.zy*B.zy;
        const real zz = A.xz*B.xz + A.yz*B.yz + A.zz*B.zz;

        return xx + yy + zz;
    }

    // Associated screw vector, operator vee
    // It assumed that the input (S) is a skew-symmetric matrix
    // Then, the output, v, is given by:
    // v = (S(3,2), S(1,3), S(2,1))

    inline __host__ __device__ real3 vee(const tensor3& S){
        return make_real3(S.zy, S.xz, S.yx);
    }

    // Pi
    // Pi(A,B) = screw(AB^T-(AB^T)^T) = screw(AB^T-BA^T)

    inline __host__ __device__ real3 Pi(const tensor3& A, const tensor3& B){
        tensor3 C = matmul(A,B.transpose());
        C = C - C.transpose();
        return vee(C);
    }




}}}
