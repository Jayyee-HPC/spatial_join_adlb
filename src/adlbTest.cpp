#define USE_UNSTABLE_GEOS_CPP_API
#include "mpi.h"
#include "adlb.h"

#include <stdio.h>
#include <stdlib.h>
#include<map>
#include<list>

#include <geos/io/WKTReader.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/index/strtree/STRtree.h>

#include "mpitype.h"
#include "MpiReadStruct.h"

#include "./spdlog/spdlog.h"
#include "./spdlog/cfg/env.h" // support for loading levels from the environment variable#include "MpiReadStruct.h"

//mpicc -o prog nq.c libadlb.a libmpigis.a -lm
//mpirun -np 2 ./prog 
void server_func(int my_world_rank)
{
	double time_t1, time_t2;
	time_t1 = MPI_Wtime();
	printf("Server Rank %d\n", my_world_rank);
	ADLB_Server(100000000000,(double)0.0);
	time_t2 = MPI_Wtime();
	printf("Server %d; Time %f.\n", my_world_rank, time_t2-time_t1);
}

void work_func(int my_world_rank, string file_path_1, string file_path_2, MPI_Comm worker_comm)
{
	if(worker_comm != MPI_COMM_NULL) 
	{
		double tv_begin = mpi_wtime();
		list<Geometry*> list_geoms_1, list_geoms_2;
		list<string> *list_strs;
		MpiReadStruct mpiReader;
		mpiReader.ReadAll(filepath1, &list_geoms_1);
		mpiReader.ReadGeomsSplit(filepath2, BLOCK_SIZE, &list_strs, worker_comm);
		mpiReader.ReadGeomsFromStr(list_strs, &list_geoms_2);
		
		spdlog::info("Rank {}, {} : {}", my_world_rank, list_geoms_1.size(), list_geoms_2->size());
		
		double tv_end_reading = mpi_wtime();

		geos::index::strtree::strtree index;

		for (list<geometry*>::iterator itr = list_geoms_1.begin() ; itr != list_geoms_1.end(); ++itr) {
			geometry* p = *itr;
			index.insert( p->getenvelopeinternal(), p );
		}

		double tv_end_indexing = mpi_wtime();
		
		for (list<geometry*>::iterator itr = list_geoms_2.begin() ; itr != list_geoms_2.end(); ++itr)
		{
			
		}
	}
	else
	{
		spdlog::error("Rank {} has an invalid worker_comm", my_world_rank);
		assert(0);
	}
}

int main(int argc, char *argv[])
{
	int my_world_rank;
	int num_world_nodes;
  
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&num_world_nodes);

	const string filepath1 = argv[2];
	const string filepath2 = argv[3];


	double start_tv;
	start_tv = MPI_Wtime(); 
	
	if (num_world_nodes < 2)
	{
		printf("** must have at least 3 ranks: one master, one worker, one server \n");
        printf("size: %ld", sizeof(Geometry));
        exit(-1);
	}
   
	MPI_Comm_rank(MPI_COMM_WORLD,&my_world_rank);
   
	spdlog::set_pattern("P" + std::to_string(my_world_rank) +  " [%H:%M:%S.%e] %v");

#ifdef DEBUG
	spdlog::set_level(spdlog::level::debug); // Set global log level to debug
#else
	spdlog::set_level(spdlog::level::info); // Set global log level to info
#endif

	MpiReadStruct mpiReader;
      
	int num_types = 1;
	int types[2] = {WORK};
	const int use_debug_server = 0;
	const int num_servers = 3;
	int am_server;
	int am_debug_server;
	MPI_Comm app_comm;
 
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	 
	int num = num_world_nodes - num_servers;
	
	if (use_debug_server) num--;
	
	int ranks[num]; 
	for(int i=0; i<num; i++) {
		ranks[i] = i;
	}
		
	MPI_Group worker_group;
	MPI_Group_incl(world_group, num, ranks, &worker_group);
        
	MPI_Comm worker_comm;
	MPI_Comm_create_group(MPI_COMM_WORLD, worker_group, 0, &worker_comm);
	
	int rc = ADLB_Init(num_servers,use_debug_server,1,num_types,types,
    &am_server,&am_debug_server,&app_comm);
				   
	if (1 != rc)
	{
		spdlog::error("ADLB_INIT_ERROR for rank {}", my_world_rank);
	}
	
	if (am_server)
	{
		server_func(my_world_rank);
		ADLB_Finalize();
		MPI_Finalize();
		exit(0);
	}
	else if (am_debug_server)
	{
		spdlog::info("Debug server on rank {}", my_world_rank);
    	ADLB_Debug_server(300.0);
    	ADLB_Finalize();
    	MPI_Finalize();
    	exit(0);
	}
	else
	{
		int num_workers = num_world_nodes - num_servers;
	
		if (use_debug_server)
			num_workers--;

		cout<<my_world_rank<<" lines  " <<wkts->size() << endl;
	
  for(list<string>::iterator itr = wkts->begin(); itr != wkts->end(); ++itr) {
			char *w = &(*itr)[0u];

			rc = adlb_put(w, itr->size(), -1, my_world_rank, work, 1);
  }
	
	delete wkts;
	wkts = nullptr;
	
	double endread_tv;
	endread_tv = mpi_wtime();

	geos::index::strtree::strtree index;

	for (list<geometry*>::iterator it = geoms.begin() ; it != geoms.end(); it++) {
		geometry* p = *it;
		index.insert( p->getenvelopeinternal(), p );
  }
	
	double endbuildindex_tv;
	endbuildindex_tv = mpi_wtime();

	int req_types[2];
	req_types[0] = work;
	req_types[1] = -1;
    
	int work_type;
	int work_prio;
	int work_handle[adlb_handle_size];
	int work_len;
	int answer_rank;
	long workcount = 0;
	
	while(1){
		geos::io::wktreader wktreader;
		list<geometry*> *g1 = new list<geometry*>;
		long startflg = 0;
		long endflg = 0;	
		rc = adlb_reserve(req_types, &work_type, &work_prio, work_handle,
                          &work_len, &answer_rank);

    if(rc == adlb_no_more_work){
      aprintf(1,"got no_more_work\n");
      break;
    }else if (rc == adlb_done_by_exhaustion){
      aprintf(1,"got exhaustion\n");
			break;
    }else if (rc < 0){
      aprintf(1,"**** reserve failed rc %d\n",rc);
      exit(-1);
    }else{
			
		}
    char *work = (char *)malloc(sizeof(char) * work_len);
		rc = adlb_get_reserved(work, work_handle);

		for(int i = 0; i < work_len; i++){
			if(work[i] != '\n'){
				endflg++;
			}else{
				string* tmpstr = new string;
				tmpstr->assign(&work[startflg], &work[endflg]);
				
				startflg = endflg + 1;
				endflg = endflg + 1;
				
				try{
					geometry* tmpgeo = null;
					tmpgeo = wktreader.read(*tmpstr).release();
				
					delete tmpstr;
					tmpstr = null;
					if(tmpgeo->isvalid())
						g1->push_back(tmpgeo);
				}catch(exception &e){
					
				}
			}			
		}

		free(work); 



		for(list<geometry*>::iterator wkritr = g1->begin(); wkritr != g1->end(); wkritr ++){
			geometry* tmpgeom = *wkritr;
						
			std::vector<void *> results;
				
			prep::preparedgeometryfactory pgf;
				
			const prep::preparedgeometry* ppdgeom = pgf.create(tmpgeom).release();
				
			index.query((*wkritr)->getenvelopeinternal(), results);

			for(vector<void *>::iterator vditr = results.begin(); vditr != results.end(); vditr ++){
				void *poly2ptr = *vditr;
				geometry* qrdgeom = (geometry*)poly2ptr;
				try{
					if(ppdgeom->intersects(qrdgeom))						
						workcount++;
				}catch(exception &e){

				}
			}

			//tmpgeom->~geometry();
			pgf.destroy(ppdgeom);

		}

    delete g1;
		g1 = null;
    if (rc == adlb_no_more_work){
      aprintf(1,"got no_more_work\n");
      break;
    }else if (rc == adlb_done_by_exhaustion){
      aprintf(1,"got exhaustion\n");
      break;
    }else if (rc < 0){
      aprintf(1,"**** reserve failed rc %d\n",rc);
      exit(-1);
     } 
  } 
        
	adlb_finalize();
    
	double end_tv;
	end_tv = mpi_wtime();
	
	
	double timeforreadfile = endread_tv - start_tv;
	double timeforbuildindex = endbuildindex_tv - endread_tv;
	double timeforgeosoperation = end_tv - endbuildindex_tv;
	
	printf("%d, %f, %f, %f \n", my_world_rank, timeforreadfile, timeforbuildindex, timeforgeosoperation);
	fflush(stdout);			
	mpi_barrier(worker_comm);
	
	int total = 0;
	
	double maxtimeforreadfile = 0;
	double maxtimeforbuldindex = 0;
	double maxtimeforgeosoperation = 0;
	double mintimeforreadfile = 0;
	double mintimeforbuldindex = 0;
	double mintimeforgeosoperation = 0;
	
	mpi_reduce(&workcount, &total, 1, mpi_int, mpi_sum, 0, worker_comm);
	mpi_reduce(&timeforreadfile, &maxtimeforreadfile, 1, mpi_double, mpi_max, 0, worker_comm);
	mpi_reduce(&timeforbuildindex, &maxtimeforbuldindex, 1, mpi_double, mpi_max, 0, worker_comm);
	mpi_reduce(&timeforgeosoperation, &maxtimeforgeosoperation, 1, mpi_double, mpi_max, 0, worker_comm);
	mpi_reduce(&timeforreadfile, &mintimeforreadfile, 1, mpi_double, mpi_min, 0, worker_comm);
	mpi_reduce(&timeforbuildindex, &mintimeforbuldindex, 1, mpi_double, mpi_min, 0, worker_comm);
	mpi_reduce(&timeforgeosoperation, &mintimeforgeosoperation, 1, mpi_double, mpi_min, 0, worker_comm);
		
	if(my_world_rank == 0){
		printf("total outputs produced = %d \n", total);
		printf("time for reading files = %f :: %f\n", maxtimeforreadfile, mintimeforreadfile);
		printf("time for building index = %f :: %f\n", maxtimeforbuldindex, mintimeforbuldindex);
		printf("time for geometry operation = %f :: %f\n", maxtimeforgeosoperation, mintimeforgeosoperation);
		fflush(stdout);
	}
		
	}


	mpi_finalize();
	return 0;
}





