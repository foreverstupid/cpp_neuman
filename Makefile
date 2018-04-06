CC=g++
CFLAGS=-Wall -O2 -g -fno-exceptions -DSHOUT
LIBS=-lm -lfftw3
LDFLAGS=

SRC=solver.cpp problem.cpp kernels.cpp vector_handler.cpp \
    string_operations.cpp
OBJ=$(SRC:%.cpp=%.o)
NAME=neuman

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.cpp $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
