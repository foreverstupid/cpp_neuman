#ifndef SOLVER_CLASS_HPP
#define SOLVER_CLASS_HPP

#include <math.h>
#include <fftw3.h>
#include "problem.hpp"
#include "vecs.h"



/* solves problem */
class Solver{
    double N;               /* first moment */
    double *C;              /* second moment */

    double *m;              /* samples of birth kernel */
    double *w;              /* samples of death kernel */
    double *mC;             /* samples of [m * C] */
    double *wC;             /* samples of [w * C] */
    double *CwC;            /* samples of [Cw * C] */
    double *w_mult_C;       /* samples od wC */

    fftw_complex *tmp_C;    /* variables for holding tmp results of */
    fftw_complex *tmp_wC;   /* convolution */
    fftw_complex *tmp_back;

    fftw_complex *fft_m;    /* fft of birth kernel */
    fftw_complex *fft_w;    /* fft of death kernel */

    fftw_plan forward_C;    /* plans for fftw3 */
    fftw_plan forward_wC;
    fftw_plan backward_mC;
    fftw_plan backward_wC;
    fftw_plan backward_CwC;

public:
    /* solve current problem */
    void solve(const Problem &problem);

    /* return solution and info about output */
    Result getResult() const
    {
        return Result(N, C);
    }

private:
    /* init new solving process */
    void init(const Problem &problem);

    /* dispose resources after solving process finish */
    void clear();

    /* get zero padded samples of birth and death kernels and for C */
    void getVectors(const Problem &problem);

    /* compute FFT of birth and death kernel */
    void getMWFFT(int n);

    /* recompute needed convolutions */
    void getConvolutions(const Problem &problem);

    /* convolving function */
    void convolve(fftw_complex *f, fftw_complex *g, fftw_plan &b,
        double *res, const Problem &problem);
};

#endif
