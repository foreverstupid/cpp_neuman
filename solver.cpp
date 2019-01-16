#include "solver.hpp"

/*====================================================================*/
/*                      ABSTRACT SOLVER METHODS                       */
/*====================================================================*/
Result AbstractSolver::solve(const Problem &p)
{
    init(p);
        /* ----------- delete this -------------- */
        //double *err = new double[p.nodes()];
        double hlp;
        //FILE *out = fopen("err_int.txt", "w");
        /* ----------- delete this -------------- */


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
            hlp = (m[j] / N - w[j] + mC[j] - p.alpha() / 2 * N *
                ((C[j] + 2) * wC[j] + CwC[j])) /
                (w[j] + p.b() - p.alpha() / 2 *
                (p.b() - p.d() - p.s() * N));
        //    err[j] = fabs(hlp - C[j]) / (fabs(C[j]) + 1e-12);
            C[j] = hlp;
        }
      //  double a = vh.getIntNorm(err, p.nodes(), p.step(), p.origin());
      //  fprintf(out, "%d %20.10lf\n", i, a);
    }

    /* correcting second moment */
    for(int i = 0; i < p.nodes(); i++){
        C[i] = (C[i] + 1);
    }

#   if defined(SHOUT) && !defined(ASCETIC)
    for(; bar_chars < BAR_WIDTH; bar_chars++){
        putchar('#');
    }
    putchar('\n');
#   endif

    clear();
   // delete[] err;
   // fclose(out);

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

    tmp_C = fftw_alloc_complex(n + 1);
    tmp_wC = fftw_alloc_complex(n + 1);
    tmp_back = fftw_alloc_complex(n + 1);
    fft_m = fftw_alloc_complex(n + 1);
    fft_w = fftw_alloc_complex(n + 1);

    forward_C = fftw_plan_dft_r2c_1d(n * 2, C, tmp_C,
        FFTW_ESTIMATE);
    forward_wC = fftw_plan_dft_r2c_1d(n * 2, w_mult_C, tmp_wC,
        FFTW_ESTIMATE);
    backward_mC = fftw_plan_dft_c2r_1d(n * 2, tmp_back, mC,
        FFTW_ESTIMATE);
    backward_wC = fftw_plan_dft_c2r_1d(n * 2, tmp_back, wC,
        FFTW_ESTIMATE);
    backward_CwC = fftw_plan_dft_c2r_1d(n * 2, tmp_back, CwC,
        FFTW_ESTIMATE);

    getMWFFT(p);
}



void SolverFFT::getMWFFT(const Problem &p)
{
    /* hold x * m(x) and x * w(x) in 3D case */
    double *tmp_m;
    double *tmp_w;

    fftw_plan m_plan;
    fftw_plan w_plan;

    if(p.dimension() == 3){
        tmp_m = new double[2 * p.nodes()];
        tmp_w = new double[2 * p.nodes()];
        double x = p.origin();

        for(int i = 0; i < 2 * p.nodes(); i++){
            tmp_m[i] = 4 * M_PI * x * m[i];
            tmp_w[i] = 4 * M_PI * x * w[i];
            x += p.step();
        }

        m_plan = fftw_plan_dft_r2c_1d(p.nodes() * 2, tmp_m, fft_m,
            FFTW_ESTIMATE);
        w_plan = fftw_plan_dft_r2c_1d(p.nodes() * 2, tmp_w, fft_w,
            FFTW_ESTIMATE);
    }else{
        m_plan = fftw_plan_dft_r2c_1d(p.nodes() * 2, m, fft_m,
            FFTW_ESTIMATE);
        w_plan = fftw_plan_dft_r2c_1d(p.nodes() * 2, w, fft_w,
            FFTW_ESTIMATE);
    }

    fftw_execute(m_plan);
    fftw_execute(w_plan);

    fftw_destroy_plan(m_plan);
    fftw_destroy_plan(w_plan);

    if(p.dimension() == 3){
        delete[] tmp_m;
        delete[] tmp_w;
    }
}



void SolverFFT::clearConvolving()
{
    fftw_free(tmp_C);
    fftw_free(tmp_wC);
    fftw_free(tmp_back);
    fftw_free(fft_m);
    fftw_free(fft_w);

    fftw_destroy_plan(forward_C);
    fftw_destroy_plan(forward_wC);
    fftw_destroy_plan(backward_mC);
    fftw_destroy_plan(backward_wC);
    fftw_destroy_plan(backward_CwC);

    fftw_cleanup();
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

    fftw_execute(forward_C);
    fftw_execute(forward_wC);

    convolve(tmp_C, fft_m, backward_mC, mC, p);
    convolve(tmp_C, fft_w, backward_wC, wC, p);
    convolve(tmp_wC, tmp_C, backward_CwC, CwC, p);

    VectorHandler::shiftLeft(mC, p.nodes(), p.nodes() / 2);
    VectorHandler::shiftLeft(wC, p.nodes(), p.nodes() / 2);
    VectorHandler::shiftLeft(CwC, p.nodes(), p.nodes() / 2);
}



void SolverFFT::convolve(const fftw_complex *f, const fftw_complex *g,
    const fftw_plan &plan, double *res, const Problem &p)
{
    double re;
    double im;

    for(int i = 0; i < p.nodes() + 1; i++){
        re = f[i][0] * g[i][0] - f[i][1] * g[i][1];
        im = f[i][0] * g[i][1] + f[i][1] * g[i][0];

        tmp_back[i][0] = re;
        tmp_back[i][1] = im;
    }

    fftw_execute(plan);

    for(int i = 0; i < p.nodes() * 2; i++){
        res[i] *= p.step() / (p.nodes() * 2);
    }
}





/*=====================================================================*/
/*                         SOLVER DHT METHODS                          */
/*=====================================================================*/
void SolverDHT::initConvolving(const Problem &p)
{
    DHTMatrix = getHankelMatrix(p.nodes(), p.step());

    Hm = new double[p.nodes()];
    Hw = new double[p.nodes()];
    HC = new double[p.nodes()];
    Hw_mult_C = new double[p.nodes()];
    tmp = new double[p.nodes()];

    VectorHandler::multiplyMatVec(DHTMatrix, m, Hm, p.nodes());
    VectorHandler::multiplyMatVec(DHTMatrix, w, Hw, p.nodes());
}



void SolverDHT::clearConvolving()
{
    delete[] DHTMatrix;
    delete[] Hm;
    delete[] Hw;
    delete[] HC;
    delete[] Hw_mult_C;
    delete[] tmp;
}



void SolverDHT::getConvolutions(const Problem &p)
{
    VectorHandler::multiplyVecs(C, w, w_mult_C, p.nodes());

    VectorHandler::multiplyMatVec(DHTMatrix, C, HC, p.nodes());
    VectorHandler::multiplyMatVec(DHTMatrix, w_mult_C, Hw_mult_C,
        p.nodes());

    convolve(HC, Hm, mC, p.nodes());
    convolve(HC, Hw, wC, p.nodes());
    convolve(HC, Hw_mult_C, CwC, p.nodes());
}



void SolverDHT::convolve(double *Hf, double *Hg, double *fg, int n)
{
    for(int i = 0; i < n; i++){
        tmp[i] = 4 * M_PI * M_PI * Hf[i] * Hg[i];
    }

    VectorHandler::multiplyMatVec(DHTMatrix, tmp, fg, n);
}



double *SolverDHT::getHankelMatrix(int n, double step)
{
    double x;
    double y = 0.0;
    double *res = new double[n * n];

    for(int i = 0; i < n; i++){
        x = 0.0;
        for(int j = 0; j < n; j++){
            res[i * n + j] = j0(x * y) * x *
                VectorHandler::weight(j, n, step);
            x += step;
        }
        y += step;
    }

    //for(int i = 0; i < n; i++){
    //    res[i] = res[i * n] = 1.0;
    //}

    return res;
}





/*
 * TEMPORARY OUT OF ORDER
 */
/*=====================================================================*/
/*                         SOLVER DHT METHODS                          */
/*=====================================================================*/
/*void SolverDHT::initConvolving(const Problem &p)
{
    Hm = new double[p.nodes()];
    Hw = new double[p.nodes()];
    HC = new double[p.nodes()];
    Hw_mult_C = new double[p.nodes()];
    tmp = new double[p.nodes()];

    cudaMalloc((void **)&kx, sizeof(double) * p.nodes());
    cudaMalloc((void **)&kb, sizeof(double) * p.nodes());

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

    cudaFree(kx);
    cudaFree(kb);
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
}*/
