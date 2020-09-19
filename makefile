FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o utils.o algo.o spmat.o group.o
	gcc -g main.o utils.o algo.o spmat.o group.o -o cluster  $(LIBS)
clean:
	rm -rf *.o cluster

main.o: main.c algo.h algo.c
	gcc -g $(FLAGS) -c main.c
spmat.o: spmat.c spmat.h group.h group.c
	gcc -g $(FLAGS) -c spmat.c
group.o: group.c group.h param.h
	gcc -g $(FLAGS) -c group.c
utils.o: utils.c utils.h spmat.h group.h spmat.c group.c
	gcc -g $(FLAGS) -c utils.c
algo.o: algo.c algo.h utils.h utils.c
	gcc -g $(FLAGS) -c algo.c