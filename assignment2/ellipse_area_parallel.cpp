#include "core/utils.h"
#include "core/cxxopts.h"
#include "core/get_time.h"
#include <iomanip>
#include <iostream>
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

void* get_points_in_ellipse_parallel(void* arg){
    if(arg==NULL){
        std::cout<<"Parameter to thread is a NULL. Terminating...\n";
        throw std::logic_error("invalid argument");
    }
    timer t;
    t.start();
    argument* as = (argument*) arg;
    float maj_radius = as->maj_radius;
    float min_radius = as->min_radius;
    uint random_seed = as->random_seed;
    unsigned long n = as->n_this_thread;
    double x_coord, y_coord;
    unsigned long ellipse_count = 0;
    for (unsigned long i = 0; i < n; i++) {
        x_coord = maj_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);

        y_coord = min_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);

        if ((sqr(x_coord)/sqr(maj_radius) + sqr(y_coord)/sqr(min_radius)) <= 1.0)
            ellipse_count++;
    }
    as->count = ellipse_count;
    as->time = t.stop();
    pthread_exit(NULL);
}

// unsigned long get_points_in_ellipse(unsigned long n, uint random_seed, float maj_radius, float min_radius) {
//   unsigned long ellipse_count = 0;
//   double x_coord, y_coord;
//   for (unsigned long i = 0; i < n; i++) {
//     x_coord = maj_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
//     y_coord = min_radius * ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
//     if ((sqr(x_coord)/sqr(maj_radius) + sqr(y_coord)/sqr(min_radius)) <= 1.0)
//       ellipse_count++;
//   }
//   return ellipse_count;
// }

void ellipse_area_calculation_serial(uint n_threads,  unsigned long n, float maj_radius, float min_radius, uint r_seed) {
  timer serial_timer;
  double time_taken = 0.0;
  uint random_seed = r_seed;
  pthread_t* threads = new pthread_t[n_threads];
  argument* args = new argument[n_threads];
  serial_timer.start();
  
  unsigned long ellipse_points = 0;
  unsigned long remainder= n%n_threads;
  for(int i=0;i<n_threads;i++){
      args[i].random_seed=r_seed+i;
      args[i].maj_radius=maj_radius;
      args[i].min_radius=min_radius;
      args[i].n_this_thread = n/n_threads;
      if(i<remainder){
          args[i].n_this_thread ++;
      }
      if(pthread_create(threads+i,NULL,get_points_in_ellipse_parallel,args+i)!=0) throw std::runtime_error("failed to start a new thread!");

  }
  std::cout << "thread_id, points_generated, ellipse_points, time_taken\n";
  for(int j=0;j<n_threads;j++){
      if(pthread_join(*(threads+j),NULL)!=0) throw std::runtime_error("failed to join a thread!");
      ellipse_points += args[j].count;
      
      std::cout << j<<", " << args[j].n_this_thread << ", "
              << args[j].count << ", " << std::setprecision(TIME_PRECISION)
              << args[j].time << "\n";
  }

  
  double area_value =
      (double)maj_radius * (double)min_radius * 4.0 * (double)ellipse_points / (double)n;

  //*------------------------------------------------------------------------
  time_taken = serial_timer.stop();
  delete[] threads;
  delete[] args;
  

  std::cout << "Total points generated : " << n << "\n";
  std::cout << "Total points in ellipse : " << ellipse_points << "\n";
  std::cout << "Result : " << std::setprecision(VAL_PRECISION) << area_value
            << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Ellipse_area_calculation",
                           "Calculate area of an ellipse using serial and parallel execution");
  options.add_options(
      "custom",
      {     
          {"nThreads", "Number of threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
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
  uint n_threads = cl_options["nThreads"].as<uint>();
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float maj_radius = cl_options["majorRad"].as<float>();
  float min_radius = cl_options["minorRad"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();

  if (n_threads < 1) {
	std::cout << "Number of threads should be equal or greater than 1. Terminating..." << std::endl;
	return 1;
  }
  if (n_points < 1) {
	std::cout << "Number of points should be equal or greater than 1. Terminating..." << std::endl;
	return 1;
  }
  if(n_threads>n_points){
      std::cout<<"Number of points should be greater than number of threads. Terminating...\n";
      return 1;
  }
  if (maj_radius <=0 || min_radius<=0) {
	std::cout << "Radius should be positive. Terminating..." << std::endl;
	return 1;
  }
  std::cout << "Number of points : " << n_points << std::endl;
  std::cout<<"Number of threads : "<<n_threads<<std::endl;
  std::cout << "Major Radius : " << maj_radius << std::endl << "Minor Radius : " << min_radius << std::endl;
  std::cout << "Random Seed : " << r_seed << std::endl;

  ellipse_area_calculation_serial(n_threads, n_points, maj_radius, min_radius, r_seed);
  return 0;
}
