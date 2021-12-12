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
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Point.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/CoordinateArraySequenceFactory.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>

#include "mpitype.h"
#include "MpiReadStruct.h"

#include "./spdlog/spdlog.h"
#include "./spdlog/cfg/env.h" // support for loading levels from the environment variable#include "MpiReadStruct.h"

#define DEBUG
#define WORK 1000
//mpicc -o prog nq.c libadlb.a libmpigis.a -lm
//mpirun -np 2 ./prog 
using namespace std;
using namespace geos::geom;
int num_types = 1;
int types[2] = {WORK};

#define TYPE_ID_UNKNOWN 0
#define TYPE_ID_POINT 1
#define TYPE_ID_LINESTRING 2
#define TYPE_ID_POLYGON 3

int get_geom_type_id(Geometry *geom)
{
	if (geom->getGeometryTypeId() == GEOS_POINT)
		return TYPE_ID_POINT;
	else if (geom->getGeometryTypeId() == GEOS_LINESTRING)
		return TYPE_ID_LINESTRING;
	else if (geom->getGeometryTypeId() == GEOS_POLYGON)
		return TYPE_ID_POLYGON;
	else
		return TYPE_ID_UNKNOWN;
}

int join_adlb_task(double* data, int data_len)
{
	if (data == nullptr || data_len <= 0)return 0;
	int num_geoms = (int)data[0];
	vector<Geometry*> geoms;
	geos::geom::CoordinateArraySequenceFactory csf;
	int curr_pos = 1 + num_geoms * 2;
	geos::geom::GeometryFactory::Ptr gf = geos::geom::GeometryFactory::create();

	for (int i = 0; i < num_geoms; ++i)
	{
		int type = (int)data[1 + i * 2];
		int size = (int)data[1 + i * 2 + 1];

		geos::geom::Geometry *geo = NULL;
        geos::geom::CoordinateArraySequence *coords_arr = new geos::geom::CoordinateArraySequence(size);

		//spdlog::debug("{} {} {}", curr_pos, size, data_len);

		for (int j = 0; j < size; ++j)
		{
			coords_arr->setAt(geos::geom::Coordinate(data[curr_pos + 2 * j], data[curr_pos + 2 * j + 1]), j);
		}
		curr_pos += (2 * size);

		if (type == TYPE_ID_POLYGON)
        {
            //POLYGON
            try
            {
                // Create non-empty LinearRing instance
                geos::geom::LinearRing ring(coords_arr, gf->getDefaultInstance());
                // Exterior (clone is required here because Polygon takes ownership)
                geo = ring.clone().release();

                geos::geom::LinearRing *exterior = dynamic_cast<geos::geom::LinearRing *>(geo);
                std::unique_ptr<geos::geom::Polygon> poly(gf->getDefaultInstance()->createPolygon(exterior, nullptr));

                //geos::geom::LinearRing* ring = new geos::geom::LinearRing(coords_arr, gf);
                //std::unique_ptr<geos::geom::Polygon> poly(gf->createPolygon(ring, nullptr));
                geo = dynamic_cast<geos::geom::Geometry *>(poly.release());

                //delete coords_arr;
            }
            catch (std::exception &e)
            {
                spdlog::error("{}", e.what());
            }
        }
        else if (type = TYPE_ID_LINESTRING)
        {
            //LINESTRING
            try
            {
                std::unique_ptr<geos::geom::LineString> line(gf->getDefaultInstance()->createLineString(coords_arr));

                geo = dynamic_cast<geos::geom::Geometry *>(line.release());
            }
            catch (std::exception &e)
            {
                spdlog::error("{}", e.what());
            }
        }
        else if (type == TYPE_ID_POINT)
        {
            //POINT
            try
            {
                std::unique_ptr<geos::geom::Point> point(gf->getDefaultInstance()->createPoint(coords_arr));
                geo = dynamic_cast<geos::geom::Geometry *>(point.release());
            }
            catch (std::exception &e)
            {
                spdlog::error("{}", e.what());
            }
        }
        else
        {
            //Unkown type, try to convert to polygon
            spdlog::error("Parsing Unknown type: {}", type);

            try
            {
                // Create non-empty LinearRing instance
                geos::geom::LinearRing ring(coords_arr, gf->getDefaultInstance());
                // Exterior (clone is required here because Polygon takes ownership)
                geo = ring.clone().release();

                geos::geom::LinearRing *exterior = dynamic_cast<geos::geom::LinearRing *>(geo);
                std::unique_ptr<geos::geom::Polygon> poly(gf->getDefaultInstance()->createPolygon(exterior, nullptr));

                //geos::geom::LinearRing* ring = new geos::geom::LinearRing(coords_arr, gf);
                //std::unique_ptr<geos::geom::Polygon> poly(gf->createPolygon(ring, nullptr));
                geo = dynamic_cast<geos::geom::Geometry *>(poly.release());
            }
            catch (std::exception &e)
            {
                spdlog::error("{}", e.what());
            }
        }

		if (geo != nullptr)
			geoms.push_back(geo);
	}

	std::unique_ptr<geos::geom::prep::PreparedGeometry> pg =
                geos::geom::prep::PreparedGeometryFactory::prepare(geoms[0]);

	int join_result = 0;

	for (int i = 1; i < geoms.size(); ++i)
		if (pg->intersects(geoms[i]))
			++join_result;

	return join_result;
}

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
	int req_types[3], work_type, work_len;
	int work_prio, work_handle[ADLB_HANDLE_SIZE], answer_rank;

	req_types[0] = WORK;
    req_types[1] = -1;

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
				//string debug_str;

				//debug_str += to_string(curr_geom->getNumPoints());
				//debug_str.push_back(' ');
				int num_geoms = 1 + query_results.size();
				int num_points = curr_geom->getNumPoints();

				//spdlog::debug("Size {}", query_results.size());
				for (int i = 0; i < query_results.size(); ++i)
				{
					num_points += ((Geometry*)(query_results[i]))->getNumPoints();
					//debug_str += to_string(((Geometry*)(query_results[i]))->getNumPoints());
					//debug_str.push_back(' ');
				}

				// one extra double to record how many geometries in this buffer
				double *temp_buf;
				temp_buf = (double*)malloc((1 + 2 * num_geoms + 2 * num_points) * sizeof(double));

				temp_buf[0]	= (double)(num_geoms);

				int curr_pos = 1 + 2 * num_geoms;

				// write the first geometry from dataset 1
				temp_buf[1] = (double)(get_geom_type_id(curr_geom));
				temp_buf[2] = (double)(curr_geom->getNumPoints());
				
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
					temp_buf[3 + 2 * i] = get_geom_type_id((Geometry*)query_results[i]);
					temp_buf[3 + 2 * i + 1] = temp_coords.get()->size();

					for (int j = 0; j < temp_coords.get()->size(); ++j)
					{
						temp_buf[curr_pos + j * 2] = temp_coords.get()->getAt(j).x;
						temp_buf[curr_pos + j * 2 + 1] = temp_coords.get()->getAt(j).y;
					}
					curr_pos += (2 * temp_coords.get()->size()); 
				}

				//spdlog::debug("PUT {} : {}", (1 + 2 * num_geoms + 2 * num_points), debug_str);
				ADLB_Put(temp_buf, (1 + 2 * num_geoms + 2 * num_points) * sizeof(double), -1, my_world_rank, WORK, 1);
			}	
		}
		
		// Workers begin to ask for tasks
		int join_result = 0;
		while (1)
		{
			int rc = ADLB_Reserve(req_types,&work_type,&work_prio,work_handle,
                          &work_len,&answer_rank);;

			// get no more work
			if (rc == ADLB_NO_MORE_WORK)
			{
				spdlog::info("Rank {}, end with no more work", my_world_rank);
				break;
			}
			else if (rc == ADLB_DONE_BY_EXHAUSTION)
			{
				spdlog::info("Rank {}, done by exhaustion", my_world_rank);
				break;
			}
			else if (rc < 0)
			{
				spdlog::error("Rank {}, reserve failed", my_world_rank);
				break;
			}

			if (work_type == WORK)
			{
				double * temp_buf;
				temp_buf = (double *)malloc(work_len);
				work_len /= sizeof(double);
				rc = ADLB_Get_reserved(temp_buf, work_handle );

				join_result += join_adlb_task(temp_buf, work_len);
				//spdlog::debug("Rank {}, get work of length {}. {} {} {}", my_world_rank, work_len, temp_buf[0], temp_buf[1], temp_buf[2]);
			}
			else /* work type mismatch*/
			{
				spdlog::error("Rank {}, unknown work type", my_world_rank);
				break;
			}
		}

		spdlog::info("Rank {} finished, join result {}", my_world_rank, join_result);
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

	const string file_path_1 = argv[2];
	const string file_path_2 = argv[3];


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
      
	const int use_debug_server = 0;
	const int num_servers = std::stoi(argv[1]);
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
		ADLB_Finalize();
	}

    
	MPI_Barrier(worker_comm);
	
	/*
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
	*/
	spdlog::debug("Rank {} finished", my_world_rank);
	MPI_Finalize();
	return 0;
}





