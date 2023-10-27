CC = g++
#MINIMUM C++17
CFLAGS = -O0 -Wall -std=c++17

#specify your variorum directory
VARIORUM_DIR ?= ../variorum/install
SYS_SAGE_DIR = ../sys-sage/install

include ${VARIORUM_DIR}/share/variorum/variorum_config.mk

INC_FLAGS = ${VARIORUM_INCLUDE_FLAGS}
LNK_FLAGS = ${VARIORUM_LINK_RPATH} ${VARIORUM_LIB_FLAGS}

all: make_dir aux ranker

stream:
	@$(CC) $(CFLAGS) -Iinclude -DSTREAM_ARRAY_SIZE=$(N) -DCACHE_LINE_SIZE=$(S) -o stream$(D) src/stream.cpp -fopenmp -Llib -lrankeraux

streamv2:
	@$(CC) -DSTREAM_ARRAY_SIZE=$(N) -o stream$(D) src/streamv2copy.cpp -fopenmp

dgemm:
	@$(CC) $(CFLAGS) $(INC_FLAGS) -Iinclude -DN=$(N) -o dgemm$(D) src/dgemm.cpp $(LNK_FLAGS) -fopenmp -Llib -lrankeraux

srandom:
	@$(CC) $(CFLAGS) -Iinclude -DSTREAM_ARRAY_SIZE=$(N) -DCACHE_LINE_SIZE=$(S) -o srandom$(D) src/srandom.cpp -fopenmp -Llib -lrankeraux

aux:
	mpic++ $(CFLAGS) $(INC_FLAGS) -Iinclude -fPIC -shared -o lib/librankeraux.so src/ranker_aux.cpp $(LNK_FLAGS)

make_dir:
	mkdir -p result
	mkdir -p result/data
	mkdir -p result/topo

ranker: lib/librankeraux.so
	mpic++ $(CFLAGS) $(INC_FLAGS) -I$(SYS_SAGE_DIR)/inc/ -Iinclude -I/usr/include/libxml2 -o ranker.o -c ranker.cpp $(LNK_FLAGS) -Llib -L$(SYS_SAGE_DIR)/lib/ -lhwloc -fopenmp -lsys-sage -lstdc++fs -lrankeraux

main: ranker.o
	mpic++ $(CFLAGS) $(INC_FLAGS) -Iinclude -I$(SYS_SAGE_DIR)/inc/ -L$(SYS_SAGE_DIR)/lib/ -o main ranker.o main.cpp $(LNK_FLAGS) -Llib -lsys-sage -lhwloc -lrankeraux