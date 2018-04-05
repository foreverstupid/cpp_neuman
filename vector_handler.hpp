#ifndef VECTOR_HANDLER_CLASS_HPP
#define VECTOR_HANDLER_CLASS_HPP

#include <fftw3.h>
#include <math.h>
#include <stdio.h>

class VectorHandler{
    int dim;

public:
    /* multiply two vectors componentwise and save result in the third
     * one */
    static void multiplyVecs(const double *f, const double *g, double *fg,
        int size);

    /* multiply square matrix and vector: b = A * x */
    static void multiplyMatVec(double *A, double *x, double *b, int n);

    /* shift vector left (not rolling) */
    static void shiftLeft(double *f, int size, int shft);

    /* store vector into the file as list of pair 'x f(x)' */
    static void storeVector(const double *f, const char *path, int size,
        double step, double origin, int accurancy);

    static double weight(int i, int n, double step)
    {
        return i == 0 || i == n - 1 ? step / 2 : step;
    }

    VectorHandler() : dim(1) {}
    VectorHandler(int dim) : dim(dim) {}

    const int getDimension() const { return dim; }

    /* get integral dot product */
    double getDot(const double *f, const double *g, int size, double step,
        double origin);

    /* return integral norm of vecor */
    double getIntNorm(const double *f, int size, double step,
        double origin);

private:
    /* TODO: make jacobian more logical and faster */
    double jacobian(double x)
    {
        return
            dim == 1 ? 1.0 :
            dim == 2 ? 2 * M_PI * x :
            4 * M_PI * x * x;
    }

    double dimCoeff()
    {
        return
            dim == 1 ? 2.0 :
            dim == 2 ? 2 * M_PI :
            4 * M_PI;
    }
};

#endif
