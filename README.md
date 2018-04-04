# Equilibrium problem solution in 1D

This project was made for the equilibrium problem solution that appears in
the Ulf Dieckmann's model of ecological systems. Neuman's method is used for
solving this problem. Second order parametrized closure is used in this
programm.

## Compiling:

### Compilation flags:

    -DSHOUT - turn off full info about iterations and turn on a progress
        bar.

    -DASCETIC - turn off all output information (even a progress bar)
        except just two numbers: the first moment and the second moment in
        zero.

    -DBAR_WIDTH=70 - set a progress bar width in characters. Default - 70.

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

    -n - set node count (positive even integer). Default: 5000.

    -r - set integration segment [-r; r] (real). Default: 20.0.

    -p - set path to store calculated data (string). Default: graph.plt.
        Note: if you want program to not create a file write "-p n".

    -a - set parameter of closure (real). Default: 0.0.

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
    s_1^m and s_1^w, second - s_2^m and s_2^w).

    e - Danchenko's exponential kerenels. Parameters: two real numbers
    (A and B parameters).

    Default kernels are kurtic.

## Example of using:

>./neuman -i 2000 -r 2.2 -kn 1.3 3.45 -p solution.txt

## Output:

    Program save calculated second moment in the file defined by the
    argument of -p in format 'x y'. Moreover it prints the first moment
    value to the standart output. If flag -DASCETIC is on, programm prints
    both the first moment value or the second moment value of zero (see
    compilation flags)

## Dependencies:

    Programm uses "fftw-3.3.7" library for convolving functions

**Author isn't responsible for mental health that can be damaged during
reading this code**
