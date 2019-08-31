#include "vector_handler.hpp"

double VectorHandler::getDot(const double *f, const double *g, int size,
    double step, double origin) const
{
    double res = 0.0;
    double x = origin;

    for(int i = 0; i < size; i++){
        res += f[i] * g[i] * weight(i, size, step) * jacobian(x);
        x += step;
    }

    return res * dimCoeff;
}



double VectorHandler::getIntNorm(const double *f, int size, double step,
    double origin) const
{
    double res = 0.0;
    double x = origin;

    for(int i = 0; i < size; i++){
        res += fabs(f[i]) * weight(i, size, step) * jacobian(x);
        x += step;
    }

    return res * dimCoeff;
}



void VectorHandler::multiplyVecs(const double *f, const double *g,
    double *fg, int size)
{
    for(int i = 0; i < size; i++){
        fg[i] = f[i] * g[i];
    }
}



void VectorHandler::multiplyMatVec(double *A, double *x, double *b, int n)
{
    double res;

    for(int i = 0; i < n; i++){
        res = 0.0;
        for(int j = 0; j < n; j++){
            res += A[i * n + j] * x[j];
        }
        b[i] = res;
    }
}



void VectorHandler::storeVector(const double *f, const char *path,
    int size, double step, double origin, int accurancy)
{
    FILE *out = fopen(path, "w");
    double x = origin;

    for(int i = 0; i < size; i++){
        fprintf(out,
            "%.*lf %.*lf\n",
            accurancy,
            x,
            accurancy,
            f[i]
        );
        x += step;
    }

    fclose(out);
}



void VectorHandler::shiftLeft(double *f, int size, int shft)
{
    for(int i = 0; i < size; i++){
        f[i] = f[i + shft];
    }
}



void VectorHandler::copy(double *dst, const double *src, int count)
{
    for(int i = 0; i < count; i++){
        dst[i] = src[i];
    }
}



double VectorHandler::getDispersion(const double *f, int size, double step,
    double origin) const
{
    double res;
    double *xs = new double[size];
    double x = origin;

    for(int i = 0; i < size; i++){
        xs[i] = x * x;
        x += step;
    }

    res = getDot(f, xs, size, step, origin);
    delete[] xs;

    return res;
}



double VectorHandler::getKurtosis(const double *f, int size, double step,
    double origin) const
{
    double res;
    double *xs = new double[size];
    double x = origin;
    double xx;
    double disp = getDispersion(f, size, step, origin);
    double norm = getIntNorm(f, size, step, origin);

    for(int i = 0; i < size; i++){
        xx = x * x;
        xs[i] = xx * xx;
        x += step;
    }

    res = getDot(f, xs, size, step, origin) / (disp * disp) * norm;
    delete[] xs;

    return res;
}
