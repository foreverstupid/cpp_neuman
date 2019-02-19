#include "problem.hpp"

const char *Problem::default_path = "graph.plt";



Problem::Problem()
{
    kernels = new KurticKernels(1.0, 1.0);

    _b = 1.0;
    _s = 1.0;
    _d = 0.0;
    _alpha = 0.0;

    _R = 20.0;
    n_count = 5000;
    i_count = 1000;
    dim = 1;

    _path = default_path;
    acc = 5;
}



int Problem::init(int argc, char **argv)
{
    int res;
    int i = 1;

    while(i < argc){
        if(argv[i][0] != '-'){
            fprintf(stderr, "### Expected argument, but get: '%s'\n",
                argv[i]);
            return expected_arg_error;
        }

        if(argv[i][1] != 'h' && !isDigit(argv[i][1]) && !argv[i + 1]){
            fprintf(stderr, "### Empty value after '%s'\n", argv[i]);
            return empty_arg_error;
        }

        if((res = handleArgument(&i, argv)) != success){
            return res;
        }
    }

    if(_R < 0.0){
        _R = getR();    
    }

    _step = _R / (n_count - 1);
    orgn = 0.0;

    return success;
}



int Problem::handleArgument(int *i, char **argv)
{
    int res;

    switch(argv[*i][1]){
        case 'n':
            n_count = str2int(argv[*i + 1]);
            break;
        case 'i':
            i_count = str2int(argv[*i + 1]);
            break;
        case 'r':
            if(argv[*i + 1][0] == 'n' && !argv[*i + 1][1]){
                _R = -1.0;
            }else{
                _R = str2double(argv[*i + 1]);
            }
            break;
        case 'D':
            dim = str2int(argv[*i + 1]);
            if(dim > 3 || dim < 1){
                fprintf(stderr, "### Wrong dimension\n");
                return dim_error;
            }
            break;
        case 'p':
            _path = argv[*i + 1];
            if(_path[0] == 'n' && !_path[1]){
                _path = 0;
            }
            break;
        case 'k':
            if((res = setKernels(i, argv)) != success){
                return res;
            }
            break;
        case 'a':
            _alpha = str2double(argv[*i + 1]);
            break;
        case 'd':
            _d = str2double(argv[*i + 1]);
            break;
        case 'b':
            _b = str2double(argv[*i + 1]);
            break;
        case 's':
            _s = str2double(argv[*i + 1]);
            break;
        case 'e':
            acc = str2int(argv[*i + 1]);
            break;
        case 'h':
            return help;
        default:
            fprintf(stderr, "### Unknown argument '%s'\n", argv[*i]);
            return unknown_arg_error;
    }

    *i += 2;
    return success;
}



int Problem::setKernels(int *i, char **argv)
{
    if(!isNumber(argv[*i + 1]) || !isNumber(argv[*i + 2])){
        fprintf(stderr, "### Invalid kernel parameters\n");
        return kernel_params_error;
    }

    switch(argv[*i][2]){
        case 'k':
            kernels = new KurticKernels(
                str2double(argv[*i + 1]),
                str2double(argv[*i + 2])
            );
            break;
        case 'K':
            kernels = new KurticKernels(
                str2double(argv[*i + 1]),
                str2double(argv[*i + 2]),
                str2double(argv[*i + 3]),
                str2double(argv[*i + 4])
            );
            *i += 2;
            break;
        case 'n':
            kernels = new NormalKernels(
                str2double(argv[*i + 1]),
                str2double(argv[*i + 2])
            );
            break;
        case 'e':
            kernels = new ExponentKernels(
                str2double(argv[*i + 1]),
                str2double(argv[*i + 2])
            );
            break;
        default:
            fprintf(stderr, "### Unknown kernel type '%s'\n", argv[*i]);
            return kernel_type_error;
    }

    *i += 1;
    return success;
}



double Problem::getDispersionM() const
{
    double res;
    double *m = new double[n_count];
    double x = orgn;
    double nm;
    VectorHandler vh = VectorHandler(dim);

    for(int i = 0; i < n_count; i++){
        m[i] = kernels->m(x);
        x += _step;
    }

    nm = vh.getIntNorm(m, n_count, _step, orgn);
    res = vh.getDispersion(m, n_count, _step, orgn);
    delete[] m;

    return res / nm;
}



double Problem::getDispersionW() const
{
    double res;
    double *w = new double[n_count];
    double x = orgn;
    double nw;
    VectorHandler vh = VectorHandler(dim);

    for(int i = 0; i < n_count; i++){
        w[i] = kernels->w(x);
        x += _step;
    }

    nw = vh.getIntNorm(w, n_count, _step, orgn);
    res = vh.getDispersion(w, n_count, _step, orgn);
    delete[] w;

    return res / nw;
}



double Problem::getExcessM() const
{
    double res;
    double *m = new double[n_count];
    double x = orgn;
    VectorHandler vh = VectorHandler(dim);

    for(int i = 0; i < n_count; i++){
        m[i] = kernels->m(x);
        x += _step;
    }

    res = vh.getKurtosis(m, n_count, _step, orgn);
    delete[] m;

    return res - 3.0;
}



double Problem::getExcessW() const
{
    double res;
    double *w = new double[n_count];
    double x = orgn;
    VectorHandler vh = VectorHandler(dim);

    for(int i = 0; i < n_count; i++){
        w[i] = kernels->w(x);
        x += _step;
    }

    res = vh.getKurtosis(w, n_count, _step, orgn);
    delete[] w;

    return res - 3.0;
}



/*
 * Gets current space size if user didn't describe it.
 */
double Problem::getR()
{
    double x = 0.0;
    double eps = 1e-7;
    double step = 1e-4;

    while(fabs(kernels->m(x)) > eps || fabs(kernels->w(x)) > eps){
        x += step;
    }

    return x;
}
