# Equilibrium problem solution in 1D, 2D and 3D

This project was made for solving the equilibrium problem with second 
order parametrized closure that appears in
the Ulf Dieckmann's model of ecological systems. Neuman's method is used for
solving this problem.

## Compiling:

### Compilation flags:

    -DSHOUT - turn off full info about iterations and turn on a progress
        bar.

    -DASCETIC - turn off all output information (even a progress bar)
        except just two numbers: the first moment and the second moment in
        zero.

    -DBAR_WIDTH=70 - set a progress bar width in characters. Default - 70.

    -DBAR_CHAR=\'#\' - set a progress bar character. Default - '#'.

## Building:
    
>make

## Cleaning:

>make clean

## CMD arguments:

    You can give arguments to the programm by cmd. Each argument name
    begins with '-' and usually contains only one letter. Argument value
    follows its name. All arguments you don't give will be set in the
    default values.

    -i - set iteration count (positive integer). Default: 1000.

    -n - set node count (positive integer). Default: 5000.

    -r - set integration segment [0; r] (real). Default: 20.0.

    -D - set dimension of space (1, 2 or 3). Default: 1.

    -p - set path to store calculated data (string). Default: graph.plt.
        Note: if you want program to not create a file write "-p n".

    -e - set accurancy of all output info in count of digits after point
        (positive integer or zero). Default: 5.

    -a - set parameter of closure (real from [0; 1]). Default: 0.0.

    -d - set environment death coeff (real). Default: 0.0.

    -s - set individual death coeff (real). Default: 1.0.

    -b - set individual birth coeff (real). Default: 1.0.

    -h - show reference about cmd parameters

    Setting kernel looks like:

    -k* [parameters]

    where * - is a letter that defines a kernel type, [parameters] - is a
    list of kernel parameters. There are three kernel types:

    n - normal kernels. Parameters: two real numbers (sigma for death and
    sigma for birth).

    k - kurtic kernels. Parameters: two real numbers (first difines
    s_1m and s_1w, second - s_2m and s_2w).

    e - Danchenko's exponential kerenels. Parameters: two real numbers
    (A and B parameters).

    Default kernels are kurtic.

    If you write the same cmd argument twice or more times, programm will
    be used the last one.

## Example of using:

>./neuman -i 2000 -r 2.2 -kn 1.3 3.45 -p solution.txt

## Output:

    Program save calculated second moment in the file defined by the
    argument of -p in format 'x y'. Moreover it prints the first moment
    value and the second moment of zero into the standart output.

## Dependencies:

    Programm uses CUDA 9.1 for DHT function convolving in 2D case and for
    FFT convolving in 1D and 3D case (cufft).

**Author isn't responsible for mental health that can be damaged during
reading this code**
