#ifndef VECTORS_OPERATIONS_MODULE_H
#define VECTORS_OPERATIONS_MODULE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double (*Func)(double);

double *get_vector(Func f, int size, double origin, double step);

//get weight of node for quadrature integration
static inline double weight(int i, int n, double step)
{
	return i == 0 || i == n - 1 ? step / 2 : step;
}

double get_dot(const double *C, const double *W, int n_count, double step);

double get_norm(Func f, int n_count, double step, double origin);

void mul(double *v, double fact, int size);

double *get_mul_vecs(double *f, double *g, int n_count);

void mul_vecs(double *f, double *g, double *res, int n_count);

void shift_left(double *f, int size, int shft);

void mul_vec_func(double *f, Func g, int n_count, double step, double orgn);

double get_diff(const double *Cn, const double *Cn_1, int n_count);

void store_vector(const double *C, const char *path, int size, double step,
    double origin, int accurancy);

#ifdef __cplusplus
}
#endif

#endif
