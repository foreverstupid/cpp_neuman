#include <stdio.h>
#include "solver.hpp"
#include "problem.hpp"
#ifdef DEBUG
#include "kernels.hpp"
#endif

#define REFERENCE_MESSAGE "List of possible cmd arguments:\n" \
"-k*   - set kernel type, where * is one of the letters:\n" \
"    n - normal kernels\n" \
"    k - kurtic kernels\n" \
"    e - exponential Danchencko kernels\n" \
"After kernel type you must write kernel parameters:\n\n" \
"    birth and death kernel dispertion for normal kernels\n" \
"    s0 and s1 parameters for kurtic kernels\n" \
"    A and B parameters for Danchencko kernels\n\n" \
"-a    - alpha parameter of closure\n" \
"-d    - environment death parameter\n" \
"-b    - kind birth parameter\n" \
"-s    - kind death parameter\n" \
"-r    - size of area\n" \
"-i    - iteration count\n" \
"-p    - path to store data (give 'n' to not create a data file)\n" \
"-n    - grid node count\n" \
"-e    - accurancy in signs after point\n" \
"-h    - show this help\n"

#ifdef DEBUG
void showArgs(const Problem &problem)
{
    const NormalKernels *kn;
    const KurticKernels *kk;
    const ExponentKernels *ke;

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
            "Kurtic kernels: sm = %.5lf, sw = %.5lf\n",
            kk->getS0(),
            kk->getS1()
        );
    }else if((ke = dynamic_cast<const ExponentKernels *>
        (&problem.getKernels())))
    {
        printf(
            "Exponent kernels: sm = %.5lf, sw = %.5lf\n",
            ke->getA(),
            ke->getB()
        );
    }

    printf(
        "R = %10.5lf\nn_count = %d\ni_count = %d\nb = %10.5lf\n"
        "s = %10.5lf\nd = %10.5lf\nalpha = %10.5lf\naccurancy = %d\n"
        "step = %10.5lf\ndimension = %d\n",
        problem.R(),
        problem.nodes(),
        problem.iters(),
        problem.b(),
        problem.s(),
        problem.d(),
        problem.alpha(),
        problem.accurancy(),
        problem.step(),
        problem.dimension()
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
            return 1;
        }
    }

#   ifdef DEBUG
    showArgs(equation);
#   endif

    AbstractSolver *solver;
    if(equation.dimension() == 1){
        solver = new SolverFFT();
    }else{
        solver = new SolverDHT();
    }
    Result answer = solver->solve(equation);

#ifdef ASCETIC
    printf("%15.*lf %15.*lf\n", equation.accurancy(), answer.N(),
    equation.accurancy(), answer.getC0());
#else
    printf("First moment: %.*lf\nC(0) = %.*lf\n",
        equation.accurancy(), answer.N(), equation.accurancy(),
        answer.getC0());
#endif

    if(equation.path()){
        VectorHandler::storeVector(answer.C(), equation.path(),
            equation.nodes(), equation.step(), -equation.R(),
            equation.accurancy());
    }

    delete solver;

    return 0;
}
