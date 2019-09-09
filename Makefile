CC=g++
CFLAGS=-Wall -O3 -Ofast -fno-exceptions -DASCETIC
LIBS=-lm -lfftw3

SRC=solver.cpp problem.cpp kernels.cpp vector_handler.cpp \
    string_operations.cpp
OBJ=$(SRC:%.cpp=%.o)
NAME=neuman

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(NAME): main.cpp $(OBJ)
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
