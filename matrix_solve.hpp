#ifndef MATRIX_SOLVE_MODULE_HPP
#define MATRIX_SOLVE_MODULE_HPP

#include <math.h>
#include <stdlib.h>
#ifndef DASCETIC
#include <stdio.h>
#endif



/*
 * Incapsulates square matrix of linear equation system for solving
 * it by Gauss method.
 */
class Matrix{
    double **data;      /* data is stored row-major */
    int size;

public:
    Matrix(){ data = 0; }
    Matrix(int size){ create(size); }
    ~Matrix(){ clear(); }

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

    /*
     * Inits a matrix by a new size (erases all stored data).
     */
    void resize(int size){ clear(); create(size); }

private:
    Matrix(const Matrix &A);
    void operator=(const Matrix &A);

    /*
     * Disposes matrix data.
     */
    void clear();

    /*
     * Creates a new data storage of a given size.
     */
    void create(int size);
};



/*
 * Solves linear equation system with the square matrix using Gauss
 * method: Ax = f.
 * Notabene: matrix and right part are changed during solving.
 */
double* solveGauss(Matrix &A, double *f);

#endif
