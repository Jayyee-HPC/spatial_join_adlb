GCC = g++
CC = mpicxx
#CC = mpecc --mpicc='mpicxx' -mpilog
#LIBB = -L/usr/local/Cellar/geos/3.5.0/lib 
#LIBB = -L/sw/geosoft/geos/lib -lgeos
LIBB = -L/usr/local/lib -lgeos
#CFLAGS = -O0 -Wall -c -std=c++0x -I/sw/geosoft/geos/include
CFLAGS = -O0 -Wall -Wfatal-errors -Wextra -c -std=c++0x -I/opt/geos/include -I/usr/local/include -Iinclude
LIBRA = ../lib/libshp.a ../lib/librtree.a ../lib/libadlb.a
#LIBRA = ../lib/libadlb.a
DEBUG =
#CFLAGS = -Wall -c -std=c++0x -I/usr/local/Cellar/geos/3.5.0/include
LFLAGS = -Wall

#-DNDEBUG add this to disable assert macro

all: allobj

seq: io

allobj: FileSplits.o MpiReadStruct.o MPI_File_Partitioner.o config.o
	$(CC) -O2 $(LFLAGS) -o allobj FileSplits.o MpiReadStruct.o MPI_File_Partitioner.o config.o $(LIBB) $(LIBRA)
	
	
types: testMPITypes.o
	$(CC) $(LFLAGS) -o types testMPITypes.o $(LIBB) $(LIBRA)
	 	
io: FileSplits.o 
	$(CC) $(LFLAGS) -o io FileSplits.o$(LIBB) $(LIBRA)

#classes
config.o: config.cpp 
	$(CC) $(CFLAGS) config.cpp
	
FileSplits.o: FileSplits.cpp
	$(CC) $(CFLAGS) $(DEBUG) FileSplits.cpp

MpiReadStruct.o: MpiReadStruct.cpp
	$(CC) $(CFLAGS) $(DEBUG) MpiReadStruct.cpp
	
MPI_File_Partitioner.o: MPI_File_Partitioner.cpp 
	$(CC) $(CFLAGS) $(DEBUG) MPI_File_Partitioner.cpp

#main class	
adlbTest:adlbTest.o MpiReadStruct.o  FileSplits.o MPI_File_Partitioner.o config.o 
	$(CC) $(LFLAGS) $(DEBUG) -o adlbTest adlbTest.o MpiReadStruct.o MPI_File_Partitioner.o FileSplits.o config.o $(LIBB) $(LIBRA)
	
adlbTest.o: adlbTest.cpp
	$(CC) $(DEBUG) $(CFLAGS) adlbTest.cpp	
	
mpiTest: mpiTest.o MpiReadStruct.o  FileSplits.o MPI_File_Partitioner.o config.o 
	$(CC) $(LFLAGS) $(DEBUG) -o mpiTest mpiTest.o MpiReadStruct.o MPI_File_Partitioner.o FileSplits.o config.o $(LIBB) $(LIBRA)

mpiTest.o: mpiTest.cpp
	$(CC) $(DEBUG) $(CFLAGS) mpiTest.cpp
	
clean:
	rm *o
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	