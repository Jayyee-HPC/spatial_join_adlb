#define USE_UNSTABLE_GEOS_CPP_API
#include "mpi.h"
#include "adlb.h"

#include <stdio.h>
#include <stdlib.h>
#include<map>
#include<list>

#include <geos/geom/Geometry.h>
#include <geos/io/WKTReader.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/index/strtree/STRtree.h>

#include "mpitype.h"
#include "MpiReadStruct.h"

#include "./spdlog/spdlog.h"
#include "./spdlog/cfg/env.h" // support for loading levels from the environment variable#include "MpiReadStruct.h"

#define WORK 1000
//mpicc -o prog nq.c libadlb.a libmpigis.a -lm
//mpirun -np 2 ./prog 
using namespace std;
using namespace geos::geom;
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
		double tv_begin = MPI_Wtime();
		list<Geometry*> list_geoms_1, list_geoms_2;
		list<string> *list_strs;
		MpiReadStruct mpiReader;
		mpiReader.ReadAll(file_path_1, &list_geoms_1);
		mpiReader.ReadGeomsSplit(file_path_2, BLOCK_SIZE, &list_strs, worker_comm);
		mpiReader.ReadGeomsFromStr(list_strs, &list_geoms_2);
		
		spdlog::info("Rank {}, {} : {}", my_world_rank, list_geoms_1.size(), list_geoms_2.size());
		
		double tv_end_reading = MPI_Wtime();

		geos::index::strtree::STRtree index;

		for (list<Geometry*>::iterator itr = list_geoms_1.begin() ; itr != list_geoms_1.end(); ++itr) {
			Geometry* p = *itr;
			index.insert( p->getEnvelopeInternal(), p );
		}

		double tv_end_indexing = MPI_Wtime();
		
		for (list<Geometry*>::iterator itr = list_geoms_2.begin() ; itr != list_geoms_2.end(); ++itr)
		{
			std::vector<void *> query_results;
			Geometry *curr_geom = *itr;
        	index.query(curr_geom->getEnvelopeInternal(), query_results);

			if (!query_results.empty())
			{
				int num_geoms = 1 + query_results.size();
				int num_points = curr_geom->getNumPoints();

				for (int i = 0; i < query_results.size(); ++i)
				{
					num_points += ((Geometry*)(query_results[i]))->getNumPoints();
				}

				// one extra double to record how many geometries in this buffer
				double *temp_buf;
				temp_buf = (double*)malloc((1 + num_geoms + 2 * num_points) * sizeof(double));

				temp_buf[0]	= (double)(1 + num_geoms);

				int curr_pos = 1 + num_geoms;

				// write the first geometry from dataset 1
				temp_buf[1] = (double)(curr_geom->getNumPoints());
				
				std::unique_ptr< CoordinateSequence > temp_coords = curr_geom->getCoordinates();

				for (int i = 0; i < temp_coords.get()->size(); ++i)
				{
					temp_buf[curr_pos + i * 2] = temp_coords.get()->getAt(i).x;
					temp_buf[curr_pos + i * 2 + 1] = temp_coords.get()->getAt(i).y;
				}

				curr_pos += (2 * curr_geom->getNumPoints());

				for (int i = 0; i < query_results.size(); ++i)
				{
					std::unique_ptr< CoordinateSequence > temp_coords = ((Geometry*)query_results[i])->getCoordinates();
					temp_buf[i+2] = temp_coords.get()->size();

					for (int j = 0; j < temp_coords.get()->size(); ++j)
					{
						temp_buf[curr_pos + j * 2] = temp_coords.get()->getAt(j).x;
						temp_buf[curr_pos + j * 2 + 1] = temp_coords.get()->getAt(j).y;
					}
					curr_pos += (2 * temp_coords.get()->size()); 
				}

				ADLB_Put(temp_buf, (num_geoms + 2 * num_points) * sizeof(double), -1, my_world_rank, WORK, 1);
			}	
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

		work_func(my_world_rank, file_path_1, file_path_2, worker_comm);
	

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





