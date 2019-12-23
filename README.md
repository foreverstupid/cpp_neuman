# Equilibrium Problem Solver

## Description

This program was made for solving of the equilibrium problem that appears in
the Ulf Dieckmann and Richard Law's model of ecological systems. The solution of
this problem can be described as a function satisfies an integral equation.
Neuman's method is used for finding such a function. The parametrized second
order closure of spatial moments is used in this equation. See
papers to get more information about the model and the problem.

## Compiling:

There is a makefile that processes the project compilation. You can use the
following compilation flags for `g++` to change some special behaviour. These
flags just define preprocess variables:

| Flag | Description | Default value (if present) |
| -- | -- | -- |
| **-DSHOUT** | Turn off full info about iterations and turn on a progress bar |
| **-DASCETIC** | Turn off all output information (even a progress bar) except the first moment value |
| **-DBAR_WIDTH=** *<int>* | Set a progress bar width in character number | 70 |
| **-DBAR_CHAR=** *<char>* | Set a progress bar character | '#' |

### Building:
```    
make
```

### Cleaning:
```
make clean
```

## CMD arguments:

You can give the programm arguments by the command line. Each argument name
begins with '-' and usually contains only one letter. The argument value
follows its name. All the arguments, you don't specify, will be set in the
default values.

### General arguments

| Argument | Description | Type | Default value |
| -- | -- | -- | -- |
| -D | Environment dimesionality | 1, 2 or 3 | 1 |
| -r | Integration area size (segment length, round or ball radius in 1D, 2D and 3D respectively). Set it with `n` if you want the size be calculated automatically | Positive real or `n` | 20.0 |
| -i | Iteration count | Positive integer | 1000 |
| -n | Grid dnde count | Positive integer | 5000 |
| -A | Alpha parameter of the closure | Real | 1.0 |
| -B | Beta parameter of the closure | Real | 1.0 |
| -G | Gamma parameter of the closure | Real | 1.0 |
| -d | Environment death rating | Real | 0.0 |
| -s | Competition rating (`d'`) | Real | 1.0 |
| -b | Birth rating | Real | 1.0 |
| -p | Path to store calculated second moment (set it with `n` if you don't want to save results in the file) | String or `n` | "graph.plt" |
| -e | Accurancy of the output in the signs after point | Positive integer | 5 |
| -h | Show the reference about cmd parameters | None | None |

### The kernel setting

Also, you can set the kernels of the needed type. The kernel argument
has the following form:
```
-k* [parameters]
```
where `*` is a letter that defines a kernel type, `[parameters]` is a
list of the kernel parameters (real numbers). The following table describes
possible kernel types:

| Letter | Kernel description | Parameters |
| -- | -- | -- |
| n | Gausian kernels | The standart deviations for the death and birth kernels |
| k | Equal kurtic kernels | s0m = s0w and s1m = s1w (two numbers) |
| K | General kurtic kernels | s0m, s1m, s0w, and s1w |
| e | Danchenko's exponential kerenels | A and B |
| r | Roughgarden kernels | sm, gamma_m, sw, and gamma_w |
| p | Exponent polynomial kernels | am, bm, aw, and bw |

Default kernels are the equal kurtic kernels with the parameters `1.0 1.0`.

### Numerical method choosing

The program can use not only nonlinear Neuman's method for solving the equation.
You can change the method, but **only in the 1D or 3D case**. Note, that using
linear methods, the program automatically sets the arguments `-A 1 -B 0 -G 0`,
making the equation linear. The user-specified values of these arguments are
ignored. Using linear methods in 2D case is ignoring.

To set the method use the following argument:
```
-m <method name>
```
where `<method name>` can be one of the folowing words:

| Name | Description |
| -- | -- |
| neuman | Nonlinear Neuman's method using series. This value is default and can be used in all dimensions |
| lneuman | Neuman's method modification for the linear case |
| nystrom | Nystrom's method using an algebraic equation system approximation |


### Example of using:
```
./neuman -i 2000 -r 2.2 -kn 1.3 3.45 -p solution.txt -r n -D 3 -s 2.0 -A 0.3
```

## Output:

Program saves the calculated second moment in the file defined by the
`-p` argument value (if it doesn't equal to `n`) in the two columns
(for `x` and for `C(x)`). Also, it prints the first moment value in the standart
output if **-DASCETIC** is set. If **-DSHOUT** is set then the program prints
the progress bar and, after completion, the first moment value and the second
moment value in zero. If neither **-DASCETIC** nor **-DSHOUT** is set then the
program prints the value of the approximate first moment on each iteration.

## Dependencies:

Programm uses [fftw-3.3.7](http://www.fftw.org/) library for performing
convolutions using Fourier method.

# **Author isn't responsible for the mental health that can be damaged during reading this code**
