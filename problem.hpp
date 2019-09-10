#ifndef PROBLEM_CLASS_HPP
#define PROBLEM_CLASS_HPP

#include <stdio.h>
#include "kernels.hpp"
#include "string_operations.hpp"
#include "vector_handler.hpp"



/* holds all input info about current problem */
class Problem{
public:
    /*
     * Solving methods for equation resolve.
     * If one of linear methods is used then A, B and G parameters are
     * ignored and asymmetric second order closure is used making
     * the problem linear.
     */
    enum Method{
        nonlinear_neuman,   /* Neuman method for nonlinear equation */
        linear_neuman,      /* Neuman method for linear equation */
        nystrom             /* Nystrom method for linear equation */
    };

private:
    static const char *default_path;    /* default path for output file */

    Kernels *kernels;   /* birth and death parameters */

    double _b;          /* birth rate */
    double _s;          /* competition rate */
    double _d;          /* death rate */

    double _alpha;      /* closure parameters */
    double _beta;
    double _gamma;

    Method _method;     /* equation solving method */

    int i_count;        /* iteration count */
    int n_count;        /* nodes count */

    double orgn;        /* origin for integration */
    double _R;          /* integration limit */
    double _step;       /* step between nodes */
    int dim;            /* dimension of space */

    int acc;            /* output accurancy */
    const char *_path;  /* path for storing data */
                        /* if _path = 0 output file won't be created */

public:
    /* comand line arguments parsing error */
    enum{
        success, help, unknown_arg_error, empty_arg_error, dim_error,
        expected_arg_error, kernel_params_error, kernel_type_error
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
    double beta() const { return _beta; }
    double gamma() const { return _gamma; }
    Method method() const { return _method; }
    int iters() const { return i_count; }
    int nodes() const { return n_count; }
    double origin() const { return orgn; }
    double R() const { return _R; }
    double step() const { return _step; }
    int dimension() const { return dim; }
    int accurancy() const { return acc; }
    const char *path() const { return _path; }

    /* some useful properties */
    double getDispersionM() const;
    double getDispersionW() const;
    double getExcessM() const;
    double getExcessW() const;

private:
    int handleArgument(int *i, char **argv);
    int setKernels(int *i, char **argv);
    double getR();
};



/* holds info about solution */
struct Result{
    double N;
    int dim;
    int n_count;
    double *C;

    Result() {}
    Result(double N, double *C, int n, int dim)
        : N(N), dim(dim), n_count(n), C(C)
    {}

    /* returns the second moment at the origin */
    double getC0() const
    {
        return C[0];
    }
};

#endif
