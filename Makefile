CC = g++
#MINIMUM C++17
CFLAGS = -O -Wall -std=c++17

#specify your variorum directory
VARIORUM_DIR ?= /u/home/sutantoj/ba/variorum/install

include ${VARIORUM_DIR}/share/variorum/variorum_config.mk

INC_FLAGS = ${VARIORUM_INCLUDE_FLAGS}
LNK_FLAGS = ${VARIORUM_LINK_RPATH} ${VARIORUM_LIB_FLAGS}

all: make_dir ranker

stream:
	@$(CC) $(CFLAGS) -Iinclude -DSTREAM_ARRAY_SIZE=$(N) -o stream src/stream.cpp -fopenmp

dgemm:
	@$(CC) $(CFLAGS) $(INC_FLAGS) -Iinclude -DN=$(N) -o dgemm$(D) src/dgemm.cpp $(LNK_FLAGS) -fopenmp

make_dir:
	mkdir -p result

ranker:
	mpic++ $(CFLAGS) $(INC_FLAGS) -Iinclude -I/usr/include/libxml2 -I../sys-sage/install/inc/ -L../sys-sage/install/lib/ -o ranker.o -c ranker.cpp $(LNK_FLAGS) -lhwloc -fopenmp -lsys-sage -lstdc++fs

main: ranker.o
	mpic++ $(CFLAGS) -I../sys-sage/install/inc/ -L../sys-sage/install/lib/ -o main ranker.o main.cpp $(LNK_FLAGS) -lsys-sage -lhwloc