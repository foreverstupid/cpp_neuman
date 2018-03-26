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

    _path = default_path;
    acc = 3;
}



int Problem::init(int argc, char **argv)
{
    int res;
    int i = 1;

    while(i < argc){
        if(argv[i][0] != '-'){
            fprintf(stderr, "### Unknown argument: '%s'\n", argv[i]);
            return error;
        }

        if(argv[i][1] != 'h' && !argv[i + 1]){
            fprintf(stderr, "### Empty value after '%s'\n", argv[i]);
            return error;
        }

        if((res = handleArgument(&i, argv)) != success){
            return res;
        }
    }

    _step = 2 * _R / (n_count - 1);
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
            _R = str2double(argv[*i + 1]);
            break;
        case 'p':
            _path = argv[*i + 1];
            if(_path[0] == 'n'){
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
            fprintf(stderr, "### Unknown argument '%s'\n",
                argv[*i]);
            return error;
    }

    *i += 2;
    return success;
}



int Problem::setKernels(int *i, char **argv)
{
    if(!is_num_arg(argv[*i + 1]) || !is_num_arg(argv[*i + 2])){
        fprintf(stderr, "### Invalid kernel parameters\n");
        return error;
    }

    switch(argv[*i][2]){
        case 'k':
            kernels = new KurticKernels(
                str2double(argv[*i + 1]),
                str2double(argv[*i + 2])
            );
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
            return error;
    }

    *i += 1;
    return success;
}