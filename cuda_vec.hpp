#ifndef CUDA_VECTOR_HPP
#define CUDA_VECTOR_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math.h>

enum{
    grid_size = 256,
    blck_size = 256
};

void cudaMatMulVec(double *kM, double *kx, double *kb, double *f,
    double *Hf, int n);

void cudaHankel(double *kx, double *kb, double *f, double *Hf,
    double step, int n);

#endif
