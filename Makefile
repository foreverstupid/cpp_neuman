CC=g++
CUDACC=nvcc
CFLAGS=-Wall -g -fno-exceptions -DSHOUT
LIBS=-lm -lfftw3 -lcuda -lcudart -lcufft
LDFLAGS=-L/usr/local/cuda/lib64
INCDIR= -I/usr/local/cuda-9.1/include/

SRC=solver.cpp problem.cpp kernels.cpp vector_handler.cpp \
    string_operations.cpp
OBJ=$(SRC:%.cpp=%.o)
NAME=neuman

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(NAME): main.cpp $(OBJ) cuda_operations.o
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

cuda_operations.o: cuda_operations.cu cuda_operations.hpp
	$(CUDACC) $(INCDIR) -c $<


ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
