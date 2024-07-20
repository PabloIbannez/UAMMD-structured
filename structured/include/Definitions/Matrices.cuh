#pragma once

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
