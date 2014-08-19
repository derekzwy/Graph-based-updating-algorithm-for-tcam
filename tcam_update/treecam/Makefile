CPP = g++
CFLGA =  -g -Wall -Wsign-compare -Wformat -std=c++0x -O2
LFLGA = -lm

SRC = rulesutils.cc treecam.cc rtrie.cc 
OBJ=$(SRC:.cc=.o)


%.o: %.cc
	${CPP} ${CFLGA} -c $^ -o $@ 

sp: ${OBJ}
	${CPP} ${LFLGA} ${CFLGA} -o treecam ${OBJ} 

all: treecam 

clean: 
	rm -f *.o
	rm -f update 
