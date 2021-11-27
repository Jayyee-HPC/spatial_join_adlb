#include"MpiReadStruct.h"

void MpiReadStruct :: ReadGeomsSplit(const string filepath, int blockSize, list<string> **listStr, MPI_Comm comm){
	int nprocs, rank; 

	MPI_Comm_size(comm, &nprocs); 
	MPI_Comm_rank(comm, &rank); 
  
	MPI_Info myinfo;
	MPI_Info_create(&myinfo);
	MPI_Info_set(myinfo, "access_style", "read_once,sequential"); 
	MPI_Info_set(myinfo, "collective_buffering", "true"); 
	MPI_Info_set(myinfo, "romio_cb_read", "enable");
  
	MPI_File fh;
	MPI_File_open(comm, filepath.c_str(), MPI_MODE_RDONLY, myinfo, &fh);
	MPI_File_Partitioner *partitioner = new MPI_File_Partitioner();
	
	FileSplits* filesplits;
	filesplits = partitioner->partitionLayer(fh, blockSize,comm);
	
	*listStr = filesplits->getContents();
	
	MPI_File_close(&fh);
}


void MpiReadStruct :: ReadGeomsFromStr(list<string> *lstr, list<Geometry*> *lgeos){		
	geos::io::WKTReader wktreader;
	for(list<string>::iterator it = lstr->begin(); it != lstr->end(); ++it){
		string tmpStr = *it;
		long startFlg = 0;
		long endFlg = 0;	
		geos::io::WKTReader wktreader;
		
		for(unsigned int i = 0; i < tmpStr.length(); i++){
			if(tmpStr[i] != '\n'){
				endFlg++;
			}else{				
				try{
				Geometry* tmpGeo = NULL;
				tmpGeo = wktreader.read(tmpStr.substr(startFlg, endFlg-startFlg)).release();
				if(tmpGeo->isValid())
						lgeos->push_back(tmpGeo);
				
					//tmpGeo->~Geometry();
				}
				catch(exception &e){
					//throw;
				}
				
				startFlg = endFlg + 1;
				endFlg = endFlg + 1;
			}			
		}
	}
}


void MpiReadStruct :: ReadAll(const string filepath, list<Geometry*>* lgeos) {
	ifstream file(filepath.c_str());
	string str;
	list<string> lstr;
	
	while (std::getline(file, str))
    {
        lstr.push_back(str);
    }
    
    ReadGeometries(&lstr, lgeos);
}

void MpiReadStruct:: ReadGeometry(string wkt, Geometry **geom) {
	geos::io::WKTReader wktreader;
	try
  	{
    	*geom = wktreader.read(wkt).release();
    }
    catch(exception &e)
  	{
  		*geom = NULL;
  	}
}

void MpiReadStruct:: ReadGeometries(string str, list<Geometry*> *lgeos) {
	char *dup = strdup(str.c_str());
    char *token = std::strtok(dup, "\n");
    while(token != NULL){
    	Geometry *geom;
    	ReadGeometry(string(token), &geom);
    	if(!isEmptyCollection(geom))
        	lgeos->push_back(geom);
        // the call is treated as a subsequent calls to strtok:
        // the function continues from where it left in previous invocation
        token = std::strtok(NULL, "\n");
    }
    free(dup);
}

void MpiReadStruct:: ReadGeometries(list<string>* lstr, list<Geometry*>* lgeos) {
    geos::io::WKTReader wktreader;
    Geometry *geom = NULL;
	
     for(list<string>::iterator it=lstr->begin();it!=lstr->end();++it) {
     	try
  		{
    		geom = wktreader.read(*it).release();
        	if(!isEmptyCollection(geom))
        		lgeos->push_back(geom);
  		}
  		catch(exception &e)
  		{
  			//Do nothing, just ignore
    		//cout<< e.what() <<endl;
  		}
     }
}

void MpiReadStruct:: ReadGeometries(list<string>* wkts, int blockSize, vector<GeomPVec>* geV) {	
	int wkt_size = wkts->size();
	int num_blocks = wkt_size/blockSize;
	
	geos::io::WKTReader wktreader;
	Geometry *geom;
   
	for(int i=0; i<=num_blocks; i++) {
   	GeomPVec gv;
   	geV->push_back(gv);
  }
   
  int block = 0;
  int counter = 0;
   
	for(list<string>::iterator it=wkts->begin();it!=wkts->end();++it) {
     	try
  		{
    		block = counter/blockSize;
    		geom = wktreader.read(*it).release();
			if(!isEmptyCollection(geom))
        		(*geV)[block].push_back(geom);
        	counter++;
  		}
  		catch(exception &e)
  		{
  			//Do nothing, just ignore
    		//cout<< e.what() <<endl;
  		}
     }
    delete(wkts);
}

bool MpiReadStruct :: isEmptyCollection(Geometry *geom)
{
    bool isEmpty = false;

    if(geom == NULL || geom->isEmpty())
      return false;
      
    for(unsigned int i=0; i<geom->getNumGeometries(); i++) {
       const Geometry *g = geom->getGeometryN(i);
      
       if(g->isEmpty())
          isEmpty = true;
    }
    
    return isEmpty;
}
//~MpiReadStruct(){
	
//}