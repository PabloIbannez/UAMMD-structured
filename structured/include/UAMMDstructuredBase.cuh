#ifndef __STRUCTURED_BASE__
#define __STRUCTURED_BASE__

//STD
#include <map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <functional>

//Third party
//#include"ThirdParty/json.hpp"

//UAMMD BASE
#include"uammd.cuh"
#include"utils/container.h"
#include"utils/quaternion.cuh"

//Define UAMMD-structured coordinates format, and forceTorque
namespace uammd{
namespace structured{

    //Define Computables
    using Computables = uammd::Interactor::Computables;

    struct ForceTorque{
        real4 force;
        real4 torque;
    };

    struct EnergyForceTorque{
        real  energy;
        real4 force;
        real4 torque;
    };

    struct ForceTorqueMagneticField{
        real4 force;
        real4 torque;
        real4 magneticField;
    };

    struct StateTransitionProbability{
        int4 tentativeState;
        real transitionProbability;
    };

// 2D Matrices
// Row-major order
#define UAMMD_SET_2D_ROW_MAJOR(matrix, rows, cols, row, col, value) (matrix)[(row) * (cols) + (col)] = (value)
#define UAMMD_GET_2D_ROW_MAJOR(matrix, rows, cols, row, col) (matrix)[(row) * (cols) + (col)]

// Column-major order
#define UAMMD_SET_2D_COL_MAJOR(matrix, rows, cols, row, col, value) (matrix)[(col) * (rows) + (row)] = (value)
#define UAMMD_GET_2D_COL_MAJOR(matrix, rows, cols, row, col) (matrix)[(col) * (rows) + (row)]

// 3D Matrices
// Row-major order
#define UAMMD_SET_3D_ROW_MAJOR(matrix, depth, rows, cols, dep, row, col, value) (matrix)[(dep) * (rows) * (cols) + (row) * (cols) + (col)] = (value)
#define UAMMD_GET_3D_ROW_MAJOR(matrix, depth, rows, cols, dep, row, col) (matrix)[(dep) * (rows) * (cols) + (row) * (cols) + (col)]

// Column-major order
#define UAMMD_SET_3D_COL_MAJOR(matrix, depth, rows, cols, dep, row, col, value) (matrix)[(col) * (rows) * (depth) + (row) * (depth) + (dep)] = (value)
#define UAMMD_GET_3D_COL_MAJOR(matrix, depth, rows, cols, dep, row, col) (matrix)[(col) * (rows) * (depth) + (row) * (depth) + (dep)]

}}

#endif
