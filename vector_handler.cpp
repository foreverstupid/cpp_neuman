#include "vector_handler.hpp"

double VectorHandler::getDot(const double *f, const double *g, int size,
    double step, double origin)
{
    double res = 0.0;
    double x = origin;

    for(int i = 0; i < size; i++){
        res += f[i] * g[i] * weight(i, size, step) * jacobian(x);
        x += step;
    }

    return res;
}



double VectorHandler::getIntNorm(const double *f, int size, double origin,
    double step)
{
    int res = 0.0;
    int x = origin;

    for(int i = 0; i < size; i++){
        res += f[i] * weight(i, size, step) * jacobian(x);
    }

    return res;
}



void VectorHandler::multiplyVecs(const double *f, const double *g,
    double *fg, int size)
{
    for(int i = 0; i < size; i++){
        fg[i] = f[i] * g[i];
    }
}


void VectorHandler::storeVector(const double *f, const char *path,
    int size, double step, double origin, int accurancy)
{
    FILE *out = fopen(path, "w");
    double x = origin;

    for(int i = 0; i < size; i++){
        fprintf(out,
            "%15.*lf %15.*lf\n",
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
