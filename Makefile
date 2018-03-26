CPPC=g++
CC=gcc
CFLAGS=-Wall -O3
LIBS=-lm -lfftw3
LDFLAGS=

CPP_SRC=solver.cpp problem.cpp kernels.cpp
CPP_OBJ=$(CPP_SRC:%.cpp=%.o)
C_SRC=str.c vecs.c
C_OBJ=$(C_SRC:%.c=%.o)
NAME=neuman

%.o: %.cpp %.hpp
	$(CPPC) $(CFLAGS) -c $< -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.cpp $(C_OBJ) $(CPP_OBJ)
	$(CPPC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(C_SRC) $(CPP_SRC)
	$(CPPC) -MM $^ > $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
