#pragma once

namespace uammd{
namespace structured{
namespace MatrixIndexing{

    inline __host__ __device__ int rowMajor2Dindex(const int& rows, const int& cols, const int& row, const int& col){
        return row * cols + col;
    }

    inline __host__ __device__ int colMajor2Dindex(const int& rows, const int& cols, const int& row, const int& col){
        return col * rows + row;
    }

    inline __host__ __device__ int rowMajor3Dindex(const int& depth, const int& rows, const int& cols, int dep, const int& row, const int& col){
        return dep * rows * cols + row * cols + col;
    }

    inline __host__ __device__ int colMajor3Dindex(const int& depth, const int& rows, const int& cols, int dep, const int& row, const int& col){
        return col * rows * depth + row * depth + dep;
    }

    // 2D row-major coordinates conversion
    inline __host__ __device__ int2 rowMajor2Dcoordinates(const int& rows, const int& cols, const int& index) {
        int row = index / cols;
        int col = index % cols;
        return make_int2(row, col);
    }

    // 2D column-major coordinates conversion
    inline __host__ __device__ int2 colMajor2Dcoordinates(const int& rows, const int& cols, const int& index) {
        int row = index % rows;
        int col = index / rows;
        return make_int2(row, col);
    }

    // 3D row-major coordinates conversion
    inline __host__ __device__ int3 rowMajor3Dcoordinates(const int& depth, const int& rows, const int& cols, const int& index) {
        int dep = index / (rows * cols);
        int remainder = index % (rows * cols);
        int row = remainder / cols;
        int col = remainder % cols;
        return make_int3(dep, row, col);
    }

    // 3D column-major coordinates conversion
    inline __host__ __device__ int3 colMajor3Dcoordinates(const int& depth, const int& rows, const int& cols, const int& index) {
        int dep = index % depth;
        int tmp = index / depth;
        int row = tmp % rows;
        int col = tmp / rows;
        return make_int3(dep, row, col);
    }

}}}
