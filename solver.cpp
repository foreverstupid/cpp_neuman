#include "solver.hpp"

Result Solver::solve(const Problem &p)
{
    init(p);
#   if defined(SHOUT) && !defined(ASCETIC)
    printf("Progress: ");
    bar_chars = 0;
    fflush(stdout);
#   endif

#   ifdef DEBUG
    VectorHandler::storeVector(m, "m.plt", p.nodes() * 2, p.step(), 0.0,
        p.accurancy());
    VectorHandler::storeVector(w, "w.plt", p.nodes() * 2, p.step(), 0.0,
        p.accurancy());
#   endif

    for(int i = 0; i < p.iters(); i++){
        N = (p.b() - p.d()) /
            (vh.getDot(C, w, p.nodes(), p.step(), -p.R()) + p.s());
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
            p.step(), 0.0, p.accurancy());
        VectorHandler::storeVector(wC, "wC.plt", p.nodes() * 2,
            p.step(), 0.0, p.accurancy());
        VectorHandler::storeVector(mC, "mC.plt", p.nodes() * 2,
            p.step(), 0.0, p.accurancy());
        VectorHandler::storeVector(CwC, "CwC.plt", p.nodes() * 2,
            p.step(), 0.0, p.accurancy());
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



void Solver::init(const Problem &p)
{
    int n = p.nodes();
    vh = VectorHandler(p.dimension());

    w = new double[n * 2];
    m = new double[n * 2];
    if(C){
        delete[] C;
    }
    C = new double[n * 2];

    mC = new double[n * 2];
    wC = new double[n * 2];
    CwC = new double[n * 2];
    w_mult_C = new double[n * 2];

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

    getVectors(p);
    getMWFFT(n);
}



void Solver::getVectors(const Problem &p)
{
    double x = -p.R();

    for(int i = 0; i < p.nodes(); i++){
        m[i] = p.b() * p.getKernels().m(x);
        C[i] = w[i] = p.s() * p.getKernels().w(x);
        x += p.step();
    }

    for(int i = p.nodes(); i < 2 * p.nodes(); i++){
        m[i] = w[i] = w_mult_C[i] = C[i] = wC[i] = mC[i] = CwC[i] = 0.0;
    }
}



void Solver::getMWFFT(int n)
{
    fftw_plan m_plan = fftw_plan_dft_r2c_1d(n * 2, m, fft_m,
        FFTW_ESTIMATE);
    fftw_plan w_plan = fftw_plan_dft_r2c_1d(n * 2, w, fft_w,
        FFTW_ESTIMATE);

    fftw_execute(m_plan);
    fftw_execute(w_plan);

    fftw_destroy_plan(m_plan);
    fftw_destroy_plan(w_plan);
}



void Solver::clear()
{
    delete[] m;
    delete[] w;
    delete[] mC;
    delete[] wC;
    delete[] CwC;
    delete[] w_mult_C;

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



void Solver::getConvolutions(const Problem &p)
{
    VectorHandler::multiplyVecs(C, w, w_mult_C, p.nodes());

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



void Solver::convolve(const fftw_complex *f, const fftw_complex *g,
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
