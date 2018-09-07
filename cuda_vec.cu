#include "cuda_vec.hpp"

__global__
void mult_kernel(double *vec, double *mtr, double *res, const int n)
{
    int tid = threadIdx.y + blockDim.y * threadIdx.y;
    double tmp = 0.0;

    if(tid < n){
        for(int i = 0; i < n; i++){
            tmp += vec[i] * mtr[i * tid + n];
        }
        res[tid] = tmp;
    }
}



void cudaMatMulVec(double *kM, double *kx, double *kb, double *f,
    double *Hf, int n)
{
    cudaMemcpy(kx, f, sizeof(double) * n, cudaMemcpyHostToDevice);
    mult_kernel<<<grid_size, blck_size>>>(kM, kx, kb, n);
    cudaMemcpy(Hf, kb, sizeof(double) * n, cudaMemcpyDeviceToHost);
}



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
    cudaMemcpy(kx, f, sizeof(double) * n, cudaMemcpyHostToDevice);
    hank_kernel<<<grid_size, blck_size>>>(kx, kb, step, n);
    cudaMemcpy(Hf, kb, sizeof(double) * n, cudaMemcpyDeviceToHost);
}
