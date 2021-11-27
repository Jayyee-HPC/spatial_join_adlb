#ifndef __MPI_READ_STRUCT_H_INCLUDED__
#define __MPI_READ_STRUCT_H_INCLUDED__

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <list>
#include <tuple>
#include <sstream>
#include <fstream>
#include <geos/geom/Geometry.h>
#include <geos/io/WKTReader.h>
#include "Parser.h"
#include "MPI_File_Partitioner.h"

typedef vector<Envelope> EnvVec;
typedef vector<Geometry> GeomVec;
typedef vector<Geometry*> GeomPVec;
// typedef char GeomWkt[168225];

using namespace std;

using namespace geos::geom;

/* noncontiguous access with a single collective I/O function */



class MpiReadStruct
{	
	public:
	void ReadGeomsFromStr(list<string> *lstr, list<Geometry*> *lgeos);
	void ReadGeomsSplit(const string filepath, int blockSize, list<string> **listStr, MPI_Comm comm);
	void ReadAll(const string filepath, list<Geometry*>* lgeos);
	void ReadGeometry(string wkt, Geometry **geom);
	void ReadGeometries(string str, list<Geometry*> *lgeos);
	void ReadGeometries(list<string>* lstr, list<Geometry*>* lgeos);
	void ReadGeometries(list<string>* wkts, int blockSize, vector<GeomPVec>* geV);
	bool isEmptyCollection(Geometry *geom);
   
	~MpiReadStruct(){
		
	};
};

#endif
