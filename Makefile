CC=g++
CUDACC=nvcc
CFLAGS=-Wall -Ofast -fno-exceptions -DASCETIC
LIBS=-lm -lfftw3# -lcuda -lcudart
LDFLAGS=-L/opt/cuda/lib64
INCDIR= -I/opt/cuda/include/

SRC=solver.cpp problem.cpp kernels.cpp vector_handler.cpp \
    string_operations.cpp
OBJ=$(SRC:%.cpp=%.o)
NAME=neuman

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(NAME): main.cpp $(OBJ)# cuda_vec.o
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

#cuda_vec.o: cuda_vec.cu cuda_vec.hpp
#	$(CUDACC) $(INCDIR) -c $<

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
