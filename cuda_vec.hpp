#ifndef CUDA_VECTOR_HPP
#define CUDA_VECTOR_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

void getMatMulVec(double *kM, double *kx, double *kb, double *f, double *Hf,
    int n);

#endif
