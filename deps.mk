str.o: str.c str.h
vecs.o: vecs.c vecs.h
solver.o: solver.cpp solver.hpp problem.hpp kernels.hpp str.h vecs.h
problem.o: problem.cpp problem.hpp kernels.hpp str.h
kernels.o: kernels.cpp kernels.hpp
