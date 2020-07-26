#ifndef SOLVER_CLASSES_HPP
#define SOLVER_CLASSES_HPP

#include <math.h>
#include <fftw3.h>
#include "problem.hpp"
#include "vector_handler.hpp"
#include "matrix_solve.hpp"
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



/* solves abstract equilibrium problem */
class AbstractSolver{
protected:
#   if defined(SHOUT) && !defined(ASCETIC)
    int bar_chars;          /* for progress bar */
#   endif
    VectorHandler vh;       /* object for vector operations */

public:
    virtual Result solve(const Problem &p) = 0;
    virtual ~AbstractSolver() {}

protected:
    /* init new solving process */
    virtual void init(const Problem &p) = 0;

    /* dispose resources after solving process finish */
    virtual void clear() = 0;
};



/* solves nonlinear equilibrium problem */
class NonlinearSolver : public AbstractSolver{
protected:
    double N;               /* first moment */
    double *C;              /* second moment samples */

    double *m;              /* birth kernel samples */
    double *w;              /* death kernel samples */

    double *mC;             /* samples of [m * C] */
    double *wC;             /* samples of [w * C] */
    double *CwC;            /* samples of [Cw * C] */
    double *w_mult_C;       /* samples of wC */

public:
    Result solve(const Problem &p);

protected:
    void init(const Problem &p);
    void clear();

    /* init data for convolving */
    virtual void initConvolving(const Problem &p) = 0;

    /* dispose data for convolving */
    virtual void clearConvolving() = 0;

    /* get zero padded samples of birth and death kernels and for C */
    void getVectors(const Problem &p);

    /* recompute needed convolutions */
    virtual void getConvolutions(const Problem &p) = 0;
};



/*
 * solves nonlinear equilibrium problem using fast Fourier transform
 * for calculating convolutions
 */
class SolverFFT : public NonlinearSolver{
    fftw_complex *tmp_C;    /* variables for holding tmp results of */
    fftw_complex *tmp_wC;   /* convolutions */
    fftw_complex *tmp_back;

    fftw_complex *fft_m;    /* fft of birth kernel */
    fftw_complex *fft_w;    /* fft of death kernel */

    fftw_plan forward_C;    /* plans for fftw3 */
    fftw_plan forward_wC;
    fftw_plan backward_mC;
    fftw_plan backward_wC;
    fftw_plan backward_CwC;

    void initConvolving(const Problem &p);
    void clearConvolving();
    void getConvolutions(const Problem &p);

    /* compute FFT of birth and death kernel */
    void getMWFFT(const Problem &p);

    /* convolving function */
    void convolve(const fftw_complex *f, const fftw_complex *g,
        const fftw_plan &plan, double *res, const Problem &p);
};



/*
 * solves nonlinear equilibrium problem using discrete Hankel
 * transform by matrix for calculating convolutions
 */
class SolverDHT : public NonlinearSolver{
    double *DHTMatrix;      /* matrix for hankel transform */

    double *Hm;             /* hankel transform of kernels */
    double *Hw;

    double *HC;             /* tmp variables for H[C] and H[wC] */
    double *Hw_mult_C;

    double *tmp;            /* help variable */

    void initConvolving(const Problem &p);
    void clearConvolving();
    void getConvolutions(const Problem &p);

    /* get matrix of DHT */
    static double *getHankelMatrix(int n, double step);

    /* make convolving vector hankel originals */
    void convolve(double *Hf, double *Hg, double *fg, int n);
};



/*
 * solves nonlinear equilibrium problem using discrete Hankel
 * transform in naive form for calculating convolutions
 */
class SolverDHTNaive : public NonlinearSolver{
    double *Hm;             /* kernel Hankel transforms */
    double *Hw;

    double *HC;             /* tmp arrays for H[C] and H[wC] */
    double *Hw_mult_C;

    double *tmp;            /* help storage */

    double *J;              /* matrix, contains J0 values */

    void initConvolving(const Problem &p);
    void clearConvolving();
    void getConvolutions(const Problem &p);

    /* get Hankel transform */
    void getHankelTransform(const double *f, double *Hf, const Problem &p);

    /* make convolving using Hankel images of functions */
    void convolve(const double *Hf, const double *Hg, double *fg,
        const Problem &p);
};



/*  solves linear equilibrium problem using Neuman method */
class LinearSolver : public AbstractSolver{
protected:
#   ifndef ASCETIC
    int shift;              /* for progress bar */
#   endif

    fftw_complex *fft_C;    /* fft of the second moment */
    fftw_complex *fft_m;    /* fft of birth kernel */

    fftw_plan forward_C;    /* plans for fftw3 */
    fftw_plan backward_mC;

    double *m;              /* samples of birth kernel */
    double *w;              /* samples of death kernel */
    double *C;              /* samples of the second moment */
    double *mC;             /* samples of [m * C] */

public:
    Result solve(const Problem &p);

protected:
    void init(const Problem &p);
    void clear();

    /* inits convolving */
    void initConvolving(const Problem &p);

    /* convolves birth kernel with the second moment */
    void convolve(const Problem &p);

    /* solves twin equation with the given parameter */
    void solveTwin(const Problem &p, double N);
};



/*
 * linear solver that uses Nystrom method (non-iterative)
 */
class NystromSolver : public AbstractSolver{
    Matrix A;       /* matrix of linear system approximation */
    double *f;      /* right parts approximation */

    double *w;      /* death kernel samples */
    double nm;      /* birth kernel norm */
    double nw;      /* death kernel norm */

public:
    Result solve(const Problem &p);

protected:
    void init(const Problem &p);
    void clear();

    /* equilibrium operator kernel in linear case */
    double kernel(const Problem &p, double x, double y)
    {
        return
            (p.b() * p.getKernels().m(x - y) +
             p.s() * p.getKernels().m(x) * p.getKernels().w(y) / nw) /
            (p.b() + p.s() / nw * p.getKernels().w(x)) / nm;
    }

    /* right part of the equation */
    double rightPart(const Problem &p, double x)
    {
        return 1.0 - (p.b() + p.s() * p.getKernels().m(x)) /
            (p.b() + p.s() * p.getKernels().w(x));
    }
};

#endif
