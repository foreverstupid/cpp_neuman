#include "matrix_solve.hpp"

void Matrix::create(int size)
{
    this->size = size;
    data = new double *[size];

    for(int i = 0; i < size; i++){
        data[i] = new double[size];
    }
}



void Matrix::clear()
{
    if(data)
    {
        for(int i = 0; i < size; i++){
            delete[] data[i];
        }

        delete[] data;
    }
}



int Matrix::findMaxAbsElementRow(int n)
{
    double tmp;
    double max = fabs(data[n][n]);
    int res = n;

    for(int i = n + 1; i < size; i++){
        tmp = fabs(data[i][n]);
        if(tmp > max){
            res = i;
            max = tmp;
        }
    }

    return res;
}



/*
 * Makes one gauss forward step (for the n-th string).
 */
void gaussStep(Matrix &A, double *f, int n)
{
    double m;
    double w = A(n, n);

    for(int i = n + 1; i < A.getSize(); i++){
        m = A(i, n) / w;

        for(int j = n + 1; j < A.getSize(); j++){
            A(i, j) -= m * A(n, j);
        }

        f[i] -= m * f[n];
        A(i, n) = 0.0;
    }
}



/*
 * Performs reverse Gauss calculation.
 */
double *reverseGauss(const Matrix &A, double *f)
{
    double hlp;
    double *x = new double[A.getSize()];

    for(int i = A.getSize() - 1; i >= 0; i--){
        hlp = 0.0;

        for(int j = A.getSize() - 1; j > i; j--){
            hlp += A(i, j) * x[j];
        }

        x[i] = (f[i] - hlp) / A(i, i);
    }

    return x;
}



double *solveGauss(Matrix &A, double *f)
{
    int m_row;

    for(int i = 0; i < A.getSize() - 1; i++){
        m_row = A.findMaxAbsElementRow(i);

        if(m_row != i){
            A.swapRows(i, m_row);
        }
#       ifndef DASCETIC
        printf("%d/%d\n", i, A.getSize());
#       endif
        gaussStep(A, f, i);
    }

    return reverseGauss(A, f);
}
