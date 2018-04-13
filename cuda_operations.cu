#include "cuda_operations.hpp"

__global__
void hank_kernel(double *vec, double *res, double step, const int n)
{
    int tid = threadIdx.y + blockDim.y * threadIdx.y;
    double tmp = 0.0;
    double x = 0.0;
    double y = step * tid;

    if(tid < n){
        for(int i = 0; i < n; i++){
            tmp += j0f(x * y) * vec[i] * x;
            x += step;
        }
        res[tid] = tmp;
    }
}



void cudaHankel(double *kx, double *kb, double *f, double *Hf,
    double step, int n)
{
    GPU_ASSERT(cudaMemcpy(kx, f, sizeof(double) * n,
        cudaMemcpyHostToDevice));
    hank_kernel<<<grid_size, blck_size>>>(kx, kb, step, n);
    GPU_ASSERT(cudaMemcpy(Hf, kb, sizeof(double) * n,
        cudaMemcpyDeviceToHost));
}



void cudaFFTForward(cufftHandle plan, double *f, double *cuda_f,
    cufftDoubleComplex *Ff, int n)
{
    GPU_ASSERT(cudaMemcpy(cuda_f, f, sizeof(double) * n,
        cudaMemcpyHostToDevice));
    cufftExecD2Z(plan, cuda_f, Ff);
}



void cudaFFTBackward(cufftHandle plan, cufftDoubleComplex *Ff,
    double *cuda_f, double *f, int n)
{
    cufftExecZ2D(plan, Ff, cuda_f);
    GPU_ASSERT(cudaMemcpy(f, cuda_f, sizeof(double) * n,
        cudaMemcpyDeviceToHost));
}



__global__
void complex_mul_kernel(const cufftDoubleComplex *f,
    const cufftDoubleComplex *g, cufftDoubleComplex *res, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    double re = f[tid].x * g[tid].y + f[tid].y * g[tid].x;
    double im = f[tid].x * g[tid].x - f[tid].y * g[tid].y;

    res[tid].x = re;
    res[tid].y = im;
}



void cudaMultiplyComplexVecs(const cufftDoubleComplex *f,
    const cufftDoubleComplex *g, cufftDoubleComplex *res, int n)
{
    complex_mul_kernel<<<grid_size, blck_size>>>(f, g, res, n);
}
