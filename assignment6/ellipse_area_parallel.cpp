#include "core/utils.h"
#include "core/cxxopts.h"
#include "core/get_time.h"
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <pthread.h>
#include <stdlib.h>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_THREADS "1"
#define DEFAULT_NUMBER_OF_POINTS "1000000000"
#define DEFAULT_MAJOR_RADIUS "2"
#define DEFAULT_MINOR_RADIUS "1"
#define DEFAULT_RANDOM_SEED "1"

struct argument{
    unsigned long n_this_thread;
    uint random_seed;
    float maj_radius;
    float min_radius;
    unsigned long count;
    double time;
};

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;
}


int main(int argc, char *argv[]) {
  // Initialize command line arguments
  MPI_Init(NULL,NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  cxxopts::Options options("Ellipse_area_calculation",
                           "Calculate area of an ellipse using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",
           cxxopts::value<unsigned long>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
	  {"majorRad", "Major radius",
	   cxxopts::value<float>()->default_value(DEFAULT_MAJOR_RADIUS)},
          {"minorRad", "Minor radius",
           cxxopts::value<float>()->default_value(DEFAULT_MINOR_RADIUS)},
          {"rSeed", "Random Seed",
           cxxopts::value<uint>()->default_value(DEFAULT_RANDOM_SEED)}
      });
  auto cl_options = options.parse(argc, argv);
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float maj_radius = cl_options["majorRad"].as<float>();
  float min_radius = cl_options["minorRad"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();

  if (n_points < 1) {
	std::cout << "Number of points should be equal or greater than 1. Terminating..." << std::endl;
	return 1;
  }
  if (maj_radius <=0 || min_radius<=0) {
	std::cout << "Radius should be positive. Terminating..." << std::endl;
	return 1;
  }
  if(world_rank==0){
    printf("Number of processes : %d",world_size);
    printf("Number of points : %ul\n",n_points);
    printf("Major Radius : %f\n",maj_radius);
    printf("Minor Radius : %f\n",min_radius);
    printf("Random Seed : %ud\n",r_seed);
  }
  timer global_timer;
  global_timer.start();
  unsigned long min_points_per_process = n_points/world_size;
  int excess_point = n_points%world_size;
  unsigned long n=min_points_per_process;
  if(world_rank<excess_point){
    n++;
  }
  //ellipse_area_calculation_serial(n_threads, n_points, maj_radius, min_radius, r_seed);
  timer t;
  t.start();

  uint random_seed = random_seed + world_rank;
  double x_coord, y_coord;
  unsigned long ellipse_count = 0;
  for (unsigned long i = 0; i < n; i++) {
      x_coord = maj_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);

      y_coord = min_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);

      if ((sqr(x_coord)/sqr(maj_radius) + sqr(y_coord)/sqr(min_radius)) <= 1.0)
          ellipse_count++;
  }
  double time = t.stop();

  unsigned long global_count = ellipse_count;

  if(world_rank == 0){
    printf("rank, points_generated, ellipse_points, time_taken\n" );
    printf("%d, %ul, %ul, %f", world_rank, n, ellipse_count, time);
    unsigned long local_count = 0;
    for(int j=1; j<world_size; j++){
      MPI_Recv(&local_count, 1, MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      global_count += local_count;
    }

    double global_time = global_timer.stop();
    double area_value =
        (double)maj_radius * (double)min_radius * 4.0 * (double)global_count / (double)n_points;
    std::cout << "Total points generated : " << n_points << "\n";
    std::cout << "Total points in ellipse : " << global_count << "\n";
    std::cout << "Result : " << std::setprecision(VAL_PRECISION) << area_value
              << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << global_time << "\n";
  }else{
    MPI_Send(&ellipse_count, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
    printf("%d, %ul, %ul, %f", world_rank, n, ellipse_count, time);
  }
  MPI_Finalize();
  return 0;
}
