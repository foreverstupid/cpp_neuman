#ifndef SOLVER_CLASS_HPP
#define SOLVER_CLASS_HPP

#include <math.h>
#include <cufft.h>
#include "problem.hpp"
#include "vector_handler.hpp"
#include "cuda_operations.hpp"
#ifdef DEBUG
#include <stdio.h>
#endif

#if defined(SHOUT) && !defined(ASCETIC)

#ifndef BAR_WIDTH
#define BAR_WIDTH 70
#endif

#ifndef BAR_CHAR
#define BAR_CHAR '#'
#endif

#endif



/* solves problem */
class AbstractSolver{
protected:
#   if defined(SHOUT) && !defined(ASCETIC)
    int bar_chars;          /* for progress bar */
#   endif
    VectorHandler vh;       /* object for vector operations */

    double N;               /* first moment */
    double *C;              /* second moment */

    double *m;              /* samples of birth kernel */
    double *w;              /* samples of death kernel */
    double *mC;             /* samples of [m * C] */
    double *wC;             /* samples of [w * C] */
    double *CwC;            /* samples of [Cw * C] */
    double *w_mult_C;       /* samples od wC */

public:
    AbstractSolver(){ C = 0; }
    virtual ~AbstractSolver(){ if(C){ delete[] C; } }

    /* solve current problem */
    Result solve(const Problem &p);

private:
    /* init new solving process */
    void init(const Problem &p);

    /* inti data for convolving */
    virtual void initConvolving(const Problem &p) = 0;

    /* dispose resources after solving process finish */
    void clear();

    /* dispose data for convolving */
    virtual void clearConvolving() = 0;

    /* get zero padded samples of birth and death kernels and for C */
    void getVectors(const Problem &p);

    /* recompute needed convolutions */
    virtual void getConvolutions(const Problem &p) = 0;
};



class SolverFFT : public AbstractSolver{
    cufftDoubleComplex *tmp_C;        /* variables for holding tmp */
    cufftDoubleComplex *tmp_wC;       /* results of convolving */
    cufftDoubleComplex *tmp_back;
    double *cuda_tmp;

    cufftDoubleComplex *fft_m;        /* fft of birth kernel */
    cufftDoubleComplex *fft_w;        /* fft of death kernel */

    cufftHandle forward_C;            /* plans for cufft */
    cufftHandle forward_wC;
    cufftHandle backward_mC;
    cufftHandle backward_wC;
    cufftHandle backward_CwC;


    void initConvolving(const Problem &p);
    void clearConvolving();
    void getConvolutions(const Problem &p);

    /* compute FFT of birth and death kernel */
    void getMWFFT(const Problem &p);

    /* convolving function */
    void convolve(const cufftDoubleComplex *f,
        const cufftDoubleComplex *g, const cufftHandle &plan, double *res,
        const Problem &p);
};



class SolverDHT : public AbstractSolver{
    double *Hm;             /* hankel transform of kernels */
    double *Hw;

    double *HC;             /* tmp variables for H[C] and H[wC] */
    double *Hw_mult_C;

    double *tmp;            /* help variables */
    double *kx;
    double *kb;


    void initConvolving(const Problem &p);
    void clearConvolving();
    void getConvolutions(const Problem &p);

    /* get hankel transform of vector */
    void getDHT(double *f, double *Hf, double step, int n);

    /* make convolving vector hankel originals */
    void convolve(double *Hf, double *Hg, double *fg, double step, int n);
};

#endif
