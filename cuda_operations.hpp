#ifndef CUDA_VECTOR_HPP
#define CUDA_VECTOR_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <math.h>

enum{
    grid_size = 256,
    blck_size = 256
};

void cudaHankel(double *kx, double *kb, double *f, double *Hf,
    double step, int n);

void cudaFFTForward(cufftHandle plan, double *f, double *cuda_f,
    cufftDoubleComplex *Ff, int n);

void cudaFFTBackward(cufftHandle plan, cufftDoubleComplex *Hf,
    double *cuda_f, double *f, int n);

void cudaMultiplyComplexVecs(const cufftDoubleComplex *f,
    const cufftDoubleComplex *g, cufftDoubleComplex *res, int n);

#endif
