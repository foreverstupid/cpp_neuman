#include "solver.hpp"

void Solver::solve(const Problem &problem)
{
    init(problem);

#   ifdef DEBUG
    store_vector(m, "m.plt", problem.nodes() * 2, problem.step(),
        -problem.R(), problem.accurancy());
    store_vector(w, "w.plt", problem.nodes() * 2, problem.step(),
        -problem.R(), problem.accurancy());
#   endif

    for(int i = 0; i < problem.iters(); i++){
        N = (problem.b() - problem.d()) /
            (get_dot(C, w, problem.nodes(), problem.step()) + problem.s());
        getConvolutions(problem);

#       ifdef DEBUG
        store_vector(C, "C.plt", problem.nodes() * 2, problem.step(),
            -problem.R(), problem.accurancy());
        store_vector(wC, "wC.plt", problem.nodes() * 2, problem.step(),
            -problem.R(), problem.accurancy());
        store_vector(mC, "mC.plt", problem.nodes() * 2, problem.step(),
            -problem.R(), problem.accurancy());
        store_vector(CwC, "CwC.plt", problem.nodes() * 2, problem.step(),
            -problem.R(), problem.accurancy());
#       endif

        for(int j = 0; j < problem.nodes(); j++){
            C[j] = (m[j] / N - w[j] + mC[j] - problem.alpha() / 2 * N *
                ((C[j] + 2) * wC[j] + CwC[j])) /
                (w[j] + problem.b() - problem.alpha() / 2 *
                (problem.b() - problem.d() - problem.s() * N));
        }
    }

    /* correcting second moment */
    for(int i = 0; i < problem.nodes(); i++){
        C[i]++;
    }

    clear();
}



void Solver::init(const Problem &problem)
{
    int n = problem.nodes();

    w = new double[n * 2];
    m = new double[n * 2];
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

    getVectors(problem);
    getMWFFT(n);
}



void Solver::getVectors(const Problem &problem)
{
    double x = -problem.R();

    for(int i = 0; i < problem.nodes(); i++){
        m[i] = problem.getKernels().m(x);
        C[i] = w[i] = problem.getKernels().w(x);
        x += problem.step();
    }

    for(int i = problem.nodes(); i < 2 * problem.nodes(); i++){
        m[i] = w[i] = w_mult_C[i] = C[i] = wC[i] = mC[i] = CwC[i] = 0.0;
    }
}



void Solver::getMWFFT(int n)
{
    fftw_plan m_plan = fftw_plan_dft_r2c_1d(n * 2, m, fft_m, FFTW_ESTIMATE);
    fftw_plan w_plan = fftw_plan_dft_r2c_1d(n * 2, w, fft_w, FFTW_ESTIMATE);

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



void Solver::getConvolutions(const Problem &problem)
{
    mul_vecs(C, w, w_mult_C, problem.nodes());

    for(int i = problem.nodes(); i < 2 * problem.nodes(); i++){
        w_mult_C[i] = C[i] = 0;
    }

#   ifdef DEBUG
    store_vector(w_mult_C, "w_C.plt", problem.nodes() * 2,
        problem.step(), -problem.R(), problem.accurancy());
#   endif

    fftw_execute(forward_C);
    fftw_execute(forward_wC);

    convolve(tmp_C, fft_m, backward_mC, mC, problem);
    convolve(tmp_C, fft_w, backward_wC, wC, problem);
    convolve(tmp_wC, tmp_C, backward_CwC, CwC, problem);

    shift_left(mC, problem.nodes(), problem.nodes() / 2);
    shift_left(wC, problem.nodes(), problem.nodes() / 2);
    shift_left(CwC, problem.nodes(), problem.nodes() / 2);
}



void Solver::convolve(fftw_complex *f, fftw_complex *g, fftw_plan &p,
    double *res, const Problem &problem)
{
    double re;
    double im;

    for(int i = 0; i < problem.nodes() + 1; i++){
        re = f[i][0] * g[i][0] - f[i][1] * g[i][1];
        im = f[i][0] * g[i][1] + f[i][1] * g[i][0];

        tmp_back[i][0] = re;
        tmp_back[i][1] = im;
    }

    fftw_execute(p);

    for(int i = 0; i < problem.nodes() * 2; i++){
        res[i] *= problem.step() / (problem.nodes() * 2);
    }
}
