CPP = g++
CFLGA =  -g -Wall -Wsign-compare -Wformat -std=c++0x -O2
LFLGA = -lm

SRC=update.cc rulesutils.cc rtrie.cc  tcam.cc
OBJ=$(SRC:.cc=.o)


%.o: %.cc
	${CPP} ${CFLGA} -c $^ -o $@ 

sp: ${OBJ}
	${CPP} ${LFLGA} ${CFLGA} -o update ${OBJ} 

all: update

clean: 
	rm -f *.o
	rm -f update 
