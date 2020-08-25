FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o utils.o algo.o spmat.o
	gcc main.o utils.o algo.o spmat.o -o cluster  $(LIBS)
clean:
	rm -rf *.o cluster

main.o: main.c utils.h
	gcc $(FLAGS) -c main.c
spmat.o: spmat.c spmat.h
	gcc $(FLAGS) -c spmat.c
utils.o: utils.c utils.h spmat.h
	gcc $(FLAGS) -c utils.c
algo.o: algo.c algo.h spmat.h
	gcc $(FLAGS) -c algo.c