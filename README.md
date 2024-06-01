# Equilibrium Problem Solver

## Description

This program is made for solving the equilibrium problem that appears in
the Ulf Dieckmann and Richard Law's population dynamics model. The solution to
this problem can be described as a function that satisfies an integral equation.
The second order closure of the third spatial moment is used in this equation:
$$
    T^{(2)}_{\alpha\beta\gamma}(x, y) =
        \dfrac{1}{\alpha + \beta} \left(
            \alpha \dfrac{C(x)C(y)}{N} +
            \beta  \dfrac{C(x)C(y - x)}{N} +
            \gamma \dfrac{C(y)C(y - x)}{N} -
            \beta N^3
        \right).
$$

See the papers to get more information about the model and the problem:
- https://user.iiasa.ac.at/~dieckman/reprints/DieckmannLaw2000.pdf
- https://user.iiasa.ac.at/~dieckman/reprints/LawDieckmann2000a.pdf
- https://user.iiasa.ac.at/~dieckman/reprints/MurrellEtal2004.pdf

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

You can give the program command line arguments. Each argument has a name
that starts with '-' and usually contains only one letter. The argument value
follows its name. All the arguments, you don't specify, will be set in the
default values.

### General arguments

| Argument | Description | Type | Default value |
| -- | -- | -- | -- |
| -D | Environment dimesionality | 1, 2 or 3 | 1 |
| -r | Integration area size (segment length or ball radius in 1D, 2D and 3D respectively). Set it with `n` if you want the size to be calculated automatically | Positive real or `n` | 20.0 |
| -i | Iterations count | Positive integer | 1000 |
| -n | Grid nodes count | Positive integer | 5000 |
| -A | Alpha parameter of the closure | Real | 1.0 |
| -B | Beta parameter of the closure | Real | 1.0 |
| -G | Gamma parameter of the closure | Real | 1.0 |
| -d | Environmental death rating | Real | 0.0 |
| -s | Competition rating (`d'`) | Real | 1.0 |
| -b | Birth rating | Real | 1.0 |
| -p | Path to store a calculated second moment (set it with `n` if you don't want to save results in the file) | String or `n` | "graph.plt" |
| -e | Accurancy of the output in decimal places | Positive integer | 5 |
| -h | Show the reference about command line parameters | None | None |

### The kernel setting

You can also specify the dispersal and competiton kernels.
The kernel argument has the following form:
```
-k* [parameters]
```
where `*` is a letter that defines a kernel type, `[parameters]` is a
list of the kernel parameters (real numbers). The dispersal kernel parameters
are defined before the competiton ones. The following table describes
possible kernel types:

| Letter | Kernel description | Parameters |
| -- | -- | -- |
| n | Gausian kernels | The standart deviations for the dispersal and competiton kernels |
| k | Equal kurtic kernels | s0m = s0w and s1m = s1w (two numbers) |
| K | General kurtic kernels | s0m, s1m, s0w, and s1w |
| e | Danchenko's exponential kerenels | A and B |
| r | Roughgarden kernels | sm, gamma_m, sw, and gamma_w |
| p | Exponent polynomial kernels | am, bm, aw, and bw |
| c | Constant kernels | Dispersal radius and competition radius |

Default kernels are the equal kurtic kernels with the parameters `1.0 1.0`.

### Numerical method choosing

By default the program uses nonlinear Neuman's method.
You can change the method, but **only in the 1D or 3D case**. When using linear
methods the program automatically sets the arguments `-A 1 -B 0 -G 0`,
making the equation linear (the user-specified values for these arguments are
ignored in that case). Using linear methods in 2D case is ignored.

To set the method use the following argument:
```
-m <method name>
```
where `<method name>` can be one of the folowing:

| Name | Description |
| -- | -- |
| neuman | Nonlinear Neuman's method. This value is default and can be used for all dimensions |
| lneuman | Neuman's method modification for the linear case |
| nystrom | Nystrom's method that uses an algebraic equation system approximation |


### Example:
```
./neuman -i 2000 -r 2.2 -kn 1.3 3.45 -p solution.txt -r n -D 3 -s 2.0 -A 0.3
```

## Output:

The program saves a calculated second moment in the file, defined by the
`-p` argument value (if it doesn't equal to `n`), in two columns
(for `x` and for `C(x)`). It also prints the first moment value in the standart
output if **-DASCETIC** is set. If **-DSHOUT** is set the program prints
the progress bar and, after completion, the first moment value and the second
moment value in zero. If neither **-DASCETIC** nor **-DSHOUT** is set then the
program prints the value of the approximate first moment on each iteration.

## Dependencies:

Programm uses [fftw-3.3.7](http://www.fftw.org/) library for performing
convolutions using Fourier method.

# **Author isn't responsible for the mental health that can be damaged during reading this code**
