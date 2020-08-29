FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o utils.o algo.o spmat.o group.o
	gcc main.o utils.o algo.o spmat.o group.o -o cluster  $(LIBS)
clean:
	rm -rf *.o cluster

main.o: main.c algo.h
	gcc $(FLAGS) -c main.c
spmat.o: spmat.c spmat.h group.h
	gcc $(FLAGS) -c spmat.c
group.o: group.c group.h
	gcc $(FLAGS) -c group.c
utils.o: utils.c utils.h spmat.h group.h
	gcc $(FLAGS) -c utils.c
algo.o: algo.c algo.h spmat.h group.h
	gcc $(FLAGS) -c algo.c