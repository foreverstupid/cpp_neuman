#ifndef PROBLEM_CLASS_HPP
#define PROBLEM_CLASS_HPP

#include <stdio.h>
#include "kernels.hpp"
#include "str.h"



/* holds all input info about current problem */
class Problem{
    static const char *default_path;    /* default path for output file */

    Kernels *kernels;   /* birth and death parameters */

    double _b;          /* birth rate */
    double _s;          /* competition rate */
    double _d;          /* death rate */

    double _alpha;      /* closure parameter */

    int i_count;        /* iteration count */
    int n_count;        /* nodes count */
    double _R;          /* integration segment */
    double _step;       /* step between nodes */

    int acc;            /* output accurancy */
    const char *_path;  /* path for storing data */
                        /* if _path = 0 output file won't be created */

public:
    enum{
        success, help, error
    };

    Problem();
    ~Problem(){ delete kernels; }

    /* initialize problem using cmd arguments */
    int init(int argc, char **argv);

    /* getters for fields */
    const Kernels &getKernels() const { return *kernels; }
    double b() const { return _b; }
    double s() const { return _s; }
    double d() const { return _d; }
    double alpha() const { return _alpha; }
    int iters() const { return i_count; }
    int nodes() const { return n_count; }
    double R() const { return _R; }
    double step() const { return _step; }
    int accurancy() const { return acc; }
    const char *path() const { return _path; }

private:
    int handleArgument(int *i, char **argv);
    int setKernels(int *i, char **argv);
};



/* holds info about solution */
class Result{
    double _N;
    const double *_C;

public:
    Result(double N, const double *C) : _N(N) { _C = C; }

    double N() const { return _N; }
    const double *C() const { return _C; }
};

#endif
