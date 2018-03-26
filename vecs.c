#include "vecs.h"

/* make vector from scalar function */
double *get_vector(Func f, int size, double origin, double step)
{
    double *res = malloc(sizeof(double) * size);
    double x = origin;
    int i;

    for(i = 0; i < size; i++){
        res[i] = f(x);
        x += step;
    }

    return res;
}



/* get dot product of two vectors */
double get_dot(const double *C, const double *W, int n_count, double step)
{
    double res = 0;
    int i;

	for(i = 0; i < n_count; i++){
		res += W[i] * C[i] * weight(i, n_count, step);
	}

    return res;
}



/* get 1-norm of function */
double get_norm(Func f, int n_count, double step, double origin)
{
    double res = 0;
	double x = origin;
    int i;

	for(i = 0; i < n_count; i++){
		res += f(x) * weight(i, n_count, step);
		x += step;
	}

    return res;
}



/* multiply vector by factor */
void mul(double *v, double fact, int size)
{
    int i;

    for(i = 0; i < size; i++){
        v[i] *= fact;
    }
}



/* get C-norm of difference between two vectors */
double get_diff(const double *Cn, const double *Cn_1, int n_count)
{
    double curr;
    double max = 0;
    int i;

    for(i = 0; i < n_count; i++){
        curr = fabs(Cn[i] - Cn_1[i]);
        if(curr > max){
            max = curr;
        }
    }

    return max;
}



/* get componentwise multipltiplication of vectors */
double *get_mul_vecs(double *f, double *g, int n_count)
{
	double *res = malloc(sizeof(double) * n_count);
	int i;

	for(i = 0; i < n_count; i++){
		res[i] = f[i] * g[i];
	}

	return res;
}



/* multiply vectors componentwise saving result in the third one */
void mul_vecs(double *f, double *g, double *res, int n_count)
{
    int i;

    for(i = 0; i < n_count; i++){
        res[i] = f[i] * g[i];
    }
}



/* multiply vector and function */
void mul_vec_func(double *f, Func g, int n_count, double step, double orgn)
{
    int i;
    double x = orgn;

    for(i = 0; i < n_count; i++){
        f[i] *= g(x);
        x += step;
    }
}


/* store vector and its component coordinates into the file */
void store_vector(const double *C, const char *path, int size, double step,
    double origin, int accurancy)
{
    FILE *out = fopen(path, "w");
    double x = origin;
    int i;

    for(i = 0; i < size; i++){
        fprintf(out,
            "%15.*lf %15.*lf\n",
            accurancy,
            x,
            accurancy,
            C[i]
        );
        x += step;
    }

    fclose(out);
}



void shift_left(double *f, int size, int shft)
{
    int i;

    for(i = 0; i < size; i++){
        f[i] = f[i + shft];
    }
}
