#include "kernels.hpp"

double Kernels::eps = 1e-12;
double NormalKernels::gauss_coeff = 1.0 / sqrt(2 * M_PI);



/* Returns essential radius of the kernels */
double Kernels::getR() const
{
    double step = 1e-5;
    double x = 0.0;

    while(fabs(m(x)) > eps || fabs(w(x)) > eps){
        x += step;
    }

    return x;
}



double NormalKernels::getR() const
{
    double sqrt2pi = sqrt(2 * M_PI);
    double logm = log(Kernels::eps * sigma_m * sqrt2pi);
    double logw = log(Kernels::eps * sigma_w * sqrt2pi);

    if (logm > 0 || logw > 0)
    {
        return 1.0;
    }

    double xm = sqrt(-2 * sigma_m * sigma_m * logm);
    double xw = sqrt(-2 * sigma_w * sigma_w * logw);

    return fmax(xm, xw);
}
