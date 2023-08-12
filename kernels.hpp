#ifndef KERNELS_CLASS_HPP
#define KERNELS_CLASS_HPP

#include <math.h>
#include <stdio.h>



/* abstract class for incapsulate birth and death kernels */
class Kernels{
protected:
    static double eps;

public:
    virtual ~Kernels() {}

    virtual double m(double x) const = 0;
    virtual double w(double x) const = 0;
    virtual double getR() const;
};



class KurticKernels : public Kernels{
    double s0m;
    double s1m;
    double s0w;
    double s1w;

public:
    KurticKernels(double s0, double s1)
        : s0m(s0), s1m(s1), s0w(s0), s1w(s1) {}
    KurticKernels(double s0m, double s1m, double s0w, double s1w)
        : s0m(s0m), s1m(s1m), s0w(s0w), s1w(s1w) {}

    double getS0m() const { return s0m; }
    double getS1m() const { return s1m; }
    double getS0w() const { return s0w; }
    double getS1w() const { return s1w; }

    double m(double x) const
    {
        double xx = x * x;
        return exp(-0.5 * (s0m * xx + s1m * xx * xx) / (1 + xx));
    }

    double w(double x) const
    {
        double xx = x * x;
        return exp(-0.5 * (s0w * xx + s1w * xx * xx) / (1 + xx));
    }
};



class NormalKernels : public Kernels{
    static double gauss_coeff;

    double sigma_m;
    double sigma_w;

public:
    NormalKernels(double sm, double sw) : sigma_m(sm), sigma_w(sw) {}

    double getSigmaM() const { return sigma_m; }
    double getSigmaW() const { return sigma_w; }

    double m(double x) const
    {
        double arg = x / sigma_m;
        return gauss_coeff / sigma_m * exp(-0.5 * arg * arg);
    }

    double w(double x) const
    {
        double arg = x / sigma_w;
        return gauss_coeff / sigma_w * exp(-0.5 * arg * arg);
    }

    double getR() const;
};



class ExponentKernels : public Kernels{
    double A;
    double B;

public:
    ExponentKernels(double A, double B) : A(A), B(B) {}

    double getA() const { return A; }
    double getB() const { return B; }

    double m(double x) const
    {
        return exp(-2 * fabs(x));
    }

    double w(double x) const
    {
        double ex = exp(-fabs(x));
        double Axx = A * x * x;

        return (ex *
            (Axx / 3 - 16.0 / 9 * A * fabs(x) + 56.0 / 27 * A + B / 3)) /
            (1 + ex * (Axx + B));
    }
};



class RoughgardenKernels : public Kernels{
    double sm;
    double gm;
    double sw;
    double gw;

public:
    RoughgardenKernels(double s, double g)
        : sm(s), gm(g), sw(s), gw(g) {}
    RoughgardenKernels(double sm, double gm, double sw, double gw)
        : sm(sm), gm(gm), sw(sw), gw(gw) {}

    double getSM() const { return sm; }
    double getGM() const { return gm; }
    double getSW() const { return sw; }
    double getGW() const { return sw; }

    double m(double x) const
    {
        return exp(-pow(fabs(x / sm), gm));
    }

    double w(double x) const
    {
        return exp(-pow(fabs(x / sw), gw));
    }
};



class ExponentPolynomialKernels : public Kernels{
    double am;
    double bm;
    double aw;
    double bw;

public:
    ExponentPolynomialKernels(double a, double b)
        : am(a), bm(b), aw(a), bw(b) {}
    ExponentPolynomialKernels(double am, double bm, double aw, double bw)
        : am(am), bm(bm), aw(aw), bw(bw) {}

    double getAM() const { return am; }
    double getBM() const { return bm; }
    double getAW() const { return aw; }
    double getBW() const { return bw; }

    double m(double x) const
    {
        double xx = x * x;
        return exp(-am * xx - bm * xx * xx);
    }

    double w(double x) const
    {
        double xx = x * x;
        return exp(-aw * xx - bw * xx * xx);
    }
};



class ConstKernels : public Kernels{
    double rm;
    double rw;
    double mValue;
    double wValue;

public:
    ConstKernels(double rm, double rw)
        : rm(rm), rw(rw), mValue(0.5 / rm), wValue(0.5 / rw) {}

    double getRM() const { return rm; }
    double getRW() const { return rw; }

    double m(double x) const
    {
        return x <= rm ? mValue : 0;
    }

    double w(double x) const
    {
        return x <= rw ? wValue : 0;
    }

    double getR() const;
};

#endif
