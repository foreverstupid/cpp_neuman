#ifndef MATRIX_SOLVE_MODULE_HPP
#define MATRIX_SOLVE_MODULE_HPP

#include <math.h>
#include <stdlib.h>



/*
 * Incapsulates square matrix of linear equation system for solving
 * it by Gauss method.
 */
class Matrix{
    double **data;      /* data is stored row-major */
    int size;

public:
    Matrix(int size);
    ~Matrix();

    double &operator()(int i, int j){ return data[i][j]; }
    const double &operator()(int i, int j) const { return data[i][j]; }
    int getSize() const { return size; }

    /*
     * Swaps i-th and j-th row.
     */
    void swapRows(int i, int j)
    {
        double *buf = data[i];
        data[i] = data[j];
        data[j] = buf;
    }

    /*
     * Returns the id of row with the greatest absolute value in n-th
     * column among the last n rows.
     */
    int findMaxAbsElementRow(int n);
};



/*
 * Solves linear equation system with the square matrix using Gauss
 * method: f = Ax.
 * Notabene: matrix is changed during solving.
 */
double* solveGauss(Matrix A, const double *f);

#endif
