#include "solver.hpp"

/*====================================================================*/
/*                      ABSTRACT SOLVER METHODS                       */
/*====================================================================*/
Result AbstractSolver::solve(const Problem &p)
{
    init(p);
#   if defined(SHOUT) && !defined(ASCETIC)
    printf("Progress: ");
    bar_chars = 0;
    fflush(stdout);
#   endif

#   ifdef DEBUG
    VectorHandler::storeVector(m, "m.plt", p.nodes() * 2, p.step(),
        p.origin(), p.accurancy());
    VectorHandler::storeVector(w, "w.plt", p.nodes() * 2, p.step(),
        p.origin(), p.accurancy());
#   endif

    for(int i = 0; i < p.iters(); i++){
        N = (p.b() - p.d()) /
            (vh.getDot(C, w, p.nodes(), p.step(), p.origin()) + p.s());
        getConvolutions(p);

#       ifndef ASCETIC
#       ifndef SHOUT
        printf("Iteration: %d\nN: %.*lf\n", i + 1,
            p.accurancy(), N);
#       else
        while(bar_chars < (double)BAR_WIDTH * i / p.iters()){
            putchar(BAR_CHAR);
            bar_chars++;
        }
        fflush(stdout);
#       endif
#       endif

#       ifdef DEBUG
        VectorHandler::storeVector(C, "C.plt", p.nodes() * 2,
            p.step(), p.origin(), p.accurancy());
        VectorHandler::storeVector(wC, "wC.plt", p.nodes() * 2,
            p.step(), p.origin(), p.accurancy());
        VectorHandler::storeVector(mC, "mC.plt", p.nodes() * 2,
            p.step(), p.origin(), p.accurancy());
        VectorHandler::storeVector(CwC, "CwC.plt", p.nodes() * 2,
            p.step(), p.origin(), p.accurancy());
#       endif

        for(int j = 0; j < p.nodes(); j++){
            C[j] = (m[j] / N - w[j] + mC[j] - p.alpha() / 2 * N *
                ((C[j] + 2) * wC[j] + CwC[j])) /
                (w[j] + p.b() - p.alpha() / 2 *
                (p.b() - p.d() - p.s() * N));
        }
    }

    /* correcting second moment */
    for(int i = 0; i < p.nodes(); i++){
        C[i]++;
    }

#   if defined(SHOUT) && !defined(ASCETIC)
    for(; bar_chars < BAR_WIDTH; bar_chars++){
        putchar('#');
    }
    putchar('\n');
#   endif

    clear();

    return Result(N, C, p.nodes(), p.dimension());
}



void AbstractSolver::init(const Problem &p)
{
    getVectors(p);
    initConvolving(p);
}



void AbstractSolver::clear()
{
    delete[] m;
    delete[] w;
    delete[] mC;
    delete[] wC;
    delete[] CwC;
    delete[] w_mult_C;

    clearConvolving();
}



void AbstractSolver::getVectors(const Problem &p)
{
    vh = VectorHandler(p.dimension());

    w = new double[p.nodes() * 2];
    m = new double[p.nodes() * 2];
    if(C){
        delete[] C;
    }
    C = new double[p.nodes() * 2];

    mC = new double[p.nodes() * 2];
    wC = new double[p.nodes() * 2];
    CwC = new double[p.nodes() * 2];
    w_mult_C = new double[p.nodes() * 2];

    double x = p.origin();

    for(int i = 0; i < p.nodes(); i++){
        m[i] = p.getKernels().m(x);
        w[i] = p.getKernels().w(x);
        x += p.step();
    }

    double nm = vh.getIntNorm(m, p.nodes(), p.step(), p.origin());
    double nw = vh.getIntNorm(w, p.nodes(), p.step(), p.origin());

#   ifdef DEBUG
    printf("nm = %15.5lf\nnw = %15.5lf\n", nm, nw);
#   endif

    for(int i = 0; i < p.nodes(); i++){
        m[i] *= p.b() / nm;
        C[i] = w[i] = w[i] * p.s() / nw;
    }

    for(int i = p.nodes(); i < 2 * p.nodes(); i++){
        m[i] = w[i] = w_mult_C[i] = C[i] = wC[i] = mC[i] = CwC[i] = 0.0;
    }
}





/*=====================================================================*/
/*                          SOLVER FFT METHODS                         */
/*=====================================================================*/
void SolverFFT::initConvolving(const Problem &p)
{
    int n = p.nodes();

    GPU_ASSERT(cudaMalloc((void **)&tmp_C,
        sizeof(cufftDoubleComplex) * (n + 1)));
    GPU_ASSERT(cudaMalloc((void **)&tmp_wC,
        sizeof(cufftDoubleComplex) * (n + 1)));
    GPU_ASSERT(cudaMalloc((void **)&tmp_back,
        sizeof(cufftDoubleComplex) * (n + 1)));
    GPU_ASSERT(cudaMalloc((void **)&fft_m,
        sizeof(cufftDoubleComplex) * (n + 1)));
    GPU_ASSERT(cudaMalloc((void **)&fft_w,
        sizeof(cufftDoubleComplex) * (n + 1)));
    GPU_ASSERT(cudaMalloc((void **)&cuda_tmp, sizeof(double) * 2 * n));

    cufftPlan1d(&forward_C, 2 * n, CUFFT_R2C, 1);
    cufftPlan1d(&forward_wC, 2 * n, CUFFT_R2C, 1);
    cufftPlan1d(&backward_mC, 2 * n, CUFFT_C2R, 1);
    cufftPlan1d(&backward_wC, 2 * n, CUFFT_C2R, 1);
    cufftPlan1d(&backward_CwC, 2 * n, CUFFT_C2R, 1);

    getMWFFT(p);
}



void SolverFFT::getMWFFT(const Problem &p)
{
    /* hold x * m(x) and x * w(x) in 3D case */
    double *tmp_m = m;
    double *tmp_w = w;
    cufftHandle plan;

    if(p.dimension() == 3){
        tmp_m = new double[2 * p.nodes()];
        tmp_w = new double[2 * p.nodes()];
        double x = p.origin();

        for(int i = 0; i < 2 * p.nodes(); i++){
            tmp_m[i] = 4 * M_PI * x * m[i];
            tmp_w[i] = 4 * M_PI * x * w[i];
            x += p.step();
        }
    }

    cufftPlan1d(&plan, 2 * p.nodes(), CUFFT_R2C, 1);
    cudaFFTForward(plan, tmp_m, cuda_tmp, fft_m, p.nodes() * 2);
    cudaFFTForward(plan, tmp_w, cuda_tmp, fft_w, p.nodes() * 2);

    if(p.dimension() == 3){
        delete[] tmp_m;
        delete[] tmp_w;
    }

    cufftDestroy(plan);
}



void SolverFFT::clearConvolving()
{
    GPU_ASSERT(cudaFree(tmp_C));
    GPU_ASSERT(cudaFree(tmp_wC));
    GPU_ASSERT(cudaFree(tmp_back));
    GPU_ASSERT(cudaFree(fft_m));
    GPU_ASSERT(cudaFree(fft_w));
    GPU_ASSERT(cudaFree(cuda_tmp));

    cufftDestroy(forward_C);
    cufftDestroy(forward_wC);
    cufftDestroy(backward_mC);
    cufftDestroy(backward_wC);
    cufftDestroy(backward_CwC);
}



void SolverFFT::getConvolutions(const Problem &p)
{
    VectorHandler::multiplyVecs(C, w, w_mult_C, p.nodes());

    if(p.dimension() == 3){
        double x = p.origin();

        for(int i = 0; i < p.nodes(); i++){
            w_mult_C[i] *= 4 * M_PI * x;
            x += p.step();
        }
    }

    for(int i = p.nodes(); i < 2 * p.nodes(); i++){
        w_mult_C[i] = C[i] = 0;
    }

#   ifdef DEBUG
    VectorHandler::storeVector(w_mult_C, "w_C.plt", p.nodes() * 2,
        p.step(), -p.R(), p.accurancy());
#   endif

    cudaFFTForward(forward_C, C, cuda_tmp, tmp_C, 2 * p.nodes());
    cudaFFTForward(forward_wC, C, cuda_tmp, tmp_C, 2 * p.nodes());

    convolve(tmp_C, fft_m, backward_mC, mC, p);
    convolve(tmp_C, fft_w, backward_wC, wC, p);
    convolve(tmp_wC, tmp_C, backward_CwC, CwC, p);

    VectorHandler::shiftLeft(mC, p.nodes(), p.nodes() / 2);
    VectorHandler::shiftLeft(wC, p.nodes(), p.nodes() / 2);
    VectorHandler::shiftLeft(CwC, p.nodes(), p.nodes() / 2);
}



void SolverFFT::convolve(const cufftDoubleComplex *f,
    const cufftDoubleComplex *g, const cufftHandle &plan, double *res,
    const Problem &p)
{
    cudaMultiplyComplexVecs(f, g, tmp_back, p.nodes() + 1);
    cudaFFTBackward(plan, tmp_back, cuda_tmp, res, 2 * p.nodes());

    for(int i = 0; i < p.nodes() * 2; i++){
        res[i] *= p.step() / (p.nodes() * 2);
    }
}





/*=====================================================================*/
/*                         SOLVER DHT METHODS                          */
/*=====================================================================*/
void SolverDHT::initConvolving(const Problem &p)
{
    Hm = new double[p.nodes()];
    Hw = new double[p.nodes()];
    HC = new double[p.nodes()];
    Hw_mult_C = new double[p.nodes()];
    tmp = new double[p.nodes()];

    GPU_ASSERT(cudaMalloc((void **)&kx, sizeof(double) * p.nodes()));
    GPU_ASSERT(cudaMalloc((void **)&kb, sizeof(double) * p.nodes()));

    getDHT(m, Hm, p.step(), p.nodes());
    getDHT(w, Hw, p.step(), p.nodes());
}



void SolverDHT::clearConvolving()
{
    delete[] Hm;
    delete[] Hw;
    delete[] HC;
    delete[] Hw_mult_C;
    delete[] tmp;

    GPU_ASSERT(cudaFree(kx));
    GPU_ASSERT(cudaFree(kb));
}



void SolverDHT::getDHT(double *f, double *Hf, double step, int n)
{
    cudaHankel(kx, kb, f, Hf, step, n);
}



void SolverDHT::getConvolutions(const Problem &p)
{
    VectorHandler::multiplyVecs(C, w, w_mult_C, p.nodes());

    getDHT(C, HC, p.step(), p.nodes());
    getDHT(w_mult_C, Hw_mult_C, p.step(), p.nodes());

    convolve(HC, Hm, mC, p.step(), p.nodes());
    convolve(HC, Hw, wC, p.step(), p.nodes());
    convolve(HC, Hw_mult_C, CwC, p.step(), p.nodes());
}



void SolverDHT::convolve(double *Hf, double *Hg, double *fg, double step,
    int n)
{
    for(int i = 0; i < n; i++){
        tmp[i] = 4 * M_PI * M_PI * Hf[i] * Hg[i];
    }

    getDHT(tmp, fg, step, n);
}
