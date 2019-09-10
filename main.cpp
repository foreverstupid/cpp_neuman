#include <stdio.h>
#include "solver.hpp"
#include "problem.hpp"
#ifdef DEBUG
#include "kernels.hpp"
#endif

#define REFERENCE_MESSAGE "EQUILIBRIUM EQUATION SOLVER\n\n"\
"This program is used for solving an integral equation that appears\n"\
"in the Ulf Dieckmann's biological model in a one kind case. This\n"\
"equation describes space moments dynamics. The program uses general\n"\
"second order closure of the third moment:\n\n"\
"              1   C(x)C(y)    C(x)C(y-x)    C(y)C(y-x)\n"\
"    T(x, y) =---(A-------- + B---------- + G---------- - BN^3)\n"\
"             A+B     N            N             N\n\n"\
"where A, B and G are alpha, beta and gamma parameters respectively.\n"\
"The Neuman method on a discrete grid is used for solving.\n\n"\
"List of possible cmd arguments:\n"\
"-k*   - set kernel type, where * is one of the letters:\n"\
"    n - normal kernels\n"\
"    k - kurtic kernels where m(x) = w(x)\n"\
"    K - general kurtic kernels\n"\
"    e - exponential Danchencko kernels\n"\
"    r - roughgarden kernels\n"\
"    p - exponent polynomial kernels\n"\
"    After kernel type you must write kernel parameters:\n"\
"      + birth and death kernel dispertion for normal kernels\n"\
"      + s0 and s1 parameters for kurtic kernels\n"\
"      + s0m, s1m, s0w and s1w parameters for general kurtic kernels\n"\
"      + A and B parameters for Danchencko kernels\n"\
"      + sm, gamma_m, sw and gamma_w parameters for roughgarden kernels\n"\
"      + am, bm, aw and bw parameters for exponent polynomial kernels\n"\
"-A - alpha parameter of second order closure\n"\
"-B - beta parameter of second order closure\n"\
"-G - gamma parameter of second order closure\n"\
"-m - equation solving method. Can be one of the following types:\n"\
"         neuman - Neuman method for nonlinear case\n"\
"         lneuman - Neuman method for linear case (LINEAR)\n"\
"         nystrom - Nystrom method (LINEAR)\n"\
"     Note that using method marked as LINEAR leads to ignoring\n"\
"     A, B and G parameters and using the asymmetric second order\n"\
"     closure (A = 1, B = G = 0). That makes equilibrium equation\n"\
"     linear one.\n"\
"-d - environment death parameter\n"\
"-b - kind birth parameter\n"\
"-s - kind death parameter\n"\
"-r - size of area (give 'n' to use autocomputed size)\n"\
"-D - dimension of space\n"\
"-i - iteration count\n"\
"-n - grid node count\n"\
"-p - path to store data (give 'n' to not create a data file)\n"\
"-e - accurancy in signs after point\n"\
"-h - show this help\n"

#ifdef DEBUG
void showArgs(const Problem &problem)
{
    const NormalKernels *kn;
    const KurticKernels *kk;
    const ExponentKernels *ke;
    const RoughgardenKernels *kr;
    const ExponentPolynomialKernels *kep;

    printf("-------------------------------------\n");
    if((kn = dynamic_cast<const NormalKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Normal kernels: sm = %.5lf, sw = %.5lf\n",
            kn->getSigmaM(),
            kn->getSigmaW()
        );
    }else if((kk = dynamic_cast<const KurticKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Kurtic kernels: sm0 = %.5lf, sm1 = %.5lf\n"
            "                sw0 = %.5lf, sw1 = %.5lf\n",
            kk->getS0m(),
            kk->getS1m(),
            kk->getS0w(),
            kk->getS1w()
        );
    }else if((ke = dynamic_cast<const ExponentKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Exponent kernels: sm = %.5lf, sw = %.5lf\n",
            ke->getA(),
            ke->getB()
        );
    }else if((kr = dynamic_cast<const RoughgardenKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Roughgarden kernels: sm = %.5lf, gamma_m = %.5lf\n"
            "                     sw = %.5lf, gamma_w = %.5lf\n",
            kr->getSM(),
            kr->getGM(),
            kr->getSW(),
            kr->getGW()
        );
    }else if((kep = dynamic_cast<const ExponentPolynomialKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Exponent polynomial kernels: sm = %.5lf, gamma_m = %.5lf\n"
            "                             sw = %.5lf, gamma_w = %.5lf\n",
            kep->getAM(),
            kep->getBM(),
            kep->getAW(),
            kep->getBW()
        );
    }



    const char *methodName;
    if(problem.method() == Problem::nonlinear_neuman){
        methodName = "neuman nonlinear";
    }else if(problem.method() == Problem::linear_neuman){
        methodName = "neuman linear";
    }else if(problem.method() == Problem::nystrom){
        methodName = "nystrom";
    }

    printf(
        "R = %10.5lf\nn_count = %d\ni_count = %d\nb = %10.5lf\n"
        "s = %10.5lf\nd = %10.5lf\nalpha = %10.5lf\nbeta = %10.5lf\n"
        "gamma = %10.5lf\naccurancy = %d\n"
        "step = %10.5lf\ndimension = %d\n"
        "method: %s\n",
        problem.R(),
        problem.nodes(),
        problem.iters(),
        problem.b(),
        problem.s(),
        problem.d(),
        problem.alpha(),
        problem.beta(),
        problem.gamma(),
        problem.accurancy(),
        problem.step(),
        problem.dimension(),
        methodName
    );
    if(problem.path()){
        printf("path = '%s'\n", problem.path());
    }
    printf("-------------------------------------\n");
}
#endif



int main(int argc, char **argv)
{
    int status;
    Problem equation;

    if((status = equation.init(argc, argv)) != Problem::success){
        if(status == Problem::help){
            printf(REFERENCE_MESSAGE);
            return 0;
        }else{
            printf("\nRun \"%s -h\" to get reference\n", argv[0]);
            return 1;
        }
    }

#   ifdef DEBUG
    showArgs(equation);
#   endif

    AbstractSolver *solver;
    if(equation.dimension() == 1 || equation.dimension() == 3){
        if(equation.method() == Problem::linear_neuman){
            solver = new LinearSolver();
        }else if(equation.method() == Problem::nystrom){
            solver = new NystromSolver();
        }else{
            solver = new SolverFFT();
        }
    }else{
        solver = new SolverDHTNaive();
    }
    Result answer = solver->solve(equation);

#   ifdef ASCETIC
    printf(
        "%15.*lf\n",
        equation.accurancy(),
        answer.N
    );

#   else
    printf("First moment: %.*lf\nC(0) = %.*lf\n",
        equation.accurancy(), answer.N, equation.accurancy(),
        answer.getC0());
#   endif

    if(equation.path()){
        VectorHandler::storeVector(answer.C, equation.path(),
            equation.nodes(), equation.step(), equation.origin(),
            equation.accurancy());
    }

    delete solver;

    return 0;
}
