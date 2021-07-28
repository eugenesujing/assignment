#include "core/utils.h"
#include "core/cxxopts.h"
#include "core/get_time.h"
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cstdio>


#define DEFAULT_GRID_SIZE "1000"
#define DEFAULT_NUMBER_OF_THREADS "1"
#define DEFAULT_CX "1"
#define DEFAULT_CY "1"
#define DEFAULT_TIME_STEPS "1000"
#define DEFAULT_MIDDLE_TEMP "600"




class TemperatureArray {
private:
	uint size;
	uint step;
	double Cx;
	double Cy;

  void assign(double *A, uint x, uint y, double newvalue) {
    A[x*size+y] = newvalue;
  };
  double read(double *A, uint x, uint y) {
    return A[x*size+y];
  };
public:
  double *CurrArray;
  double *PrevArray;
	TemperatureArray(uint input_size, double iCx, double iCy, double init_temp) { // create array of dimension sizexsize
		size = input_size;
 		Cx = iCx;
    		Cy = iCy;
		step = 0;
		CurrArray = (double *)malloc(size*size*sizeof(double));
    if(CurrArray==NULL){
      std::cout<<"currarray malloc failed,which may be due to the memory limit\n";
       exit(1);
    }
		PrevArray = (double *)malloc(size*size*sizeof(double));
    if(PrevArray==NULL){
      std::cout<<"prevarray malloc failed, which may be due to the memory limit\n";
      free(CurrArray);
      exit(1);
    }
		for (uint i = 0; i < size; i++)
			for (uint j = 0; j < size; j++) {
				if ((i > size/3) && (i < 2*size/3) && (j > size/3) && (j < 2*size/3)) {
					assign(PrevArray, i, j, init_temp); assign(CurrArray, i, j, init_temp);
				}
				else {
					assign(PrevArray, i, j, 0); assign (CurrArray, i, j, 0);
				}
			}
	};

	~TemperatureArray() {
		free (PrevArray);   free (CurrArray);
	};

	void IncrementStepCount() { step ++; };

	uint ReadStepCount() { return(step); };

	void ComputeNewTemp(uint x, uint y) {
		if ((x > 0) && (x < size-1) && (y > 0) && (y < size-1))
			assign(CurrArray, x, y , read(PrevArray,x,y)
				+ Cx * (read(PrevArray, x-1, y) + read(PrevArray, x+1, y) - 2*read(PrevArray, x, y))
				+ Cy * (read(PrevArray, x, y-1) + read(PrevArray, x, y+1) - 2*read(PrevArray, x, y)));
	};

	void SwapArrays() {
		double *temp = PrevArray;
		PrevArray = CurrArray;
		CurrArray = temp;
	};

	double temp(uint x, uint y) {
		return read(CurrArray, x, y);
	};
};



int main(int argc, char *argv[]) {
  MPI_Init(NULL,NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Initialize command line arguments
  cxxopts::Options options("Heat_transfer_calculation",
                           "Model heat transfer in a grid using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"gSize", "Grid Size",
           cxxopts::value<uint>()->default_value(DEFAULT_GRID_SIZE)},
          {"mTemp", "Temperature in middle of array",
           cxxopts::value<double>()->default_value(DEFAULT_MIDDLE_TEMP)},
	        {"iCX", "Coefficient of horizontal heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CX)},
          {"iCY", "Coefficient of vertical heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CY)},
          {"tSteps", "Time Steps",
           cxxopts::value<uint>()->default_value(DEFAULT_TIME_STEPS)}
      });
  auto cl_options = options.parse(argc, argv);

  uint grid_size = cl_options["gSize"].as<uint>();
   if (grid_size <=0) {
	std::cout << "Grid_size should be positive. Terminating..." << "\n";
	return 1;
  }
  if(world_size>grid_size){
    std::cout<<"Grid size should be greater than number of processes. Terminating...\n";
    return 1;
  }
  double init_temp = cl_options["mTemp"].as<double>();
  double Cx = cl_options["iCX"].as<double>();
  double Cy = cl_options["iCY"].as<double>();
  uint steps = cl_options["tSteps"].as<uint>();
   if (steps < 1) {
	std::cout << "Number of steps should be equal or greater than 1. Terminating..." << "\n";
	return 1;
  }
  if(world_rank == 0){
    std::cout << "Number of threads : " << world_size <<"\n";
    std::cout << "Grid Size : " << grid_size << "x" << grid_size << "\n";

    std::cout << "Cx : " << Cx << "\n" << "Cy : " << Cy << "\n";
    std::cout << "Temperature in the middle of grid : " << init_temp << "\n";
    std::cout << "Time Steps : " << steps << "\n";

    std::cout << "Initializing Temperature Array..." << "\n";
    std::cout << "thread_id, start_column, end_column, time_taken\n";
  }

  TemperatureArray *T = new TemperatureArray(grid_size, Cx, Cy, init_temp);
  timer global_timer;
  global_timer.start();
  // std::cout<<"Finishing initializing array...\n";
  if (!T) {
      std::cout << "Cannot Initialize Temperature Array...Terminating" << "\n";
      return 2;
  }
  int min_columns = grid_size / world_size;
  int excess_columns = grid_size % world_size;
  int startx = 0;
  int endx = 0;
   if (world_rank < excess_columns) {
      startx = world_rank * (min_columns + 1);
      endx = startx + min_columns;
   }
   else {
     startx = (excess_columns * (min_columns + 1)) + ((world_rank-excess_columns) * min_columns);
     endx = startx + min_columns - 1;
   }
   int local_stepcount;
   timer t1;
   t1.start();
   for(local_stepcount = 1; local_stepcount < steps; local_stepcount++) {
    //Compute the Temperature Array values Curr[][] in the slice allocated to this thread from Prev[][]
    for (uint x = startx; x <= endx; x++) {
		  for (uint y = 0; y < grid_size; y++) {
			  T->ComputeNewTemp(x, y);
		  }
	  }
    // --- synchronization: Send and Receive boundary columns from neighbors
    // Even processes communicate with right proces first
    // Odd  processes communicate with left process first
    if (world_rank % 2 == 0)  {   // even rank
        if (world_rank < world_size - 1)  {  // not last process
            //Send my column "end" to the right process world_rank+1
            MPI_Send(T->CurrArray +(endx * grid_size), grid_size, MPI_DOUBLE, world_rank+1, endx, MPI_COMM_WORLD);
            //Receive column "end+1" from the right process world_rank+1, populate local Curr Array
            MPI_Recv(T->CurrArray + (endx+1)*grid_size, grid_size, MPI_DOUBLE, world_rank+1, endx+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (world_rank > 0)   {  // not first process
            //Receive column "start-1" from the left process world_rank-1, populate local Curr Array
            MPI_Recv(T->CurrArray +((startx-1) * grid_size), grid_size, MPI_DOUBLE, world_rank-1, startx-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //Send my column "start" to the left process world_rank-1
            MPI_Send(T->CurrArray + startx*grid_size, grid_size, MPI_DOUBLE, world_rank-1, startx, MPI_COMM_WORLD);
        }
        }
       // even rank
    else {  // odd rank
        if (world_rank > 0)   {  // not first process
            //Receive column "start-1" from the left process world_rank-1, populate local Curr Array
            MPI_Recv(T->CurrArray +((startx-1) * grid_size), grid_size, MPI_DOUBLE, world_rank-1, startx-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //Send my column "start" to the left process world_rank-1
            MPI_Send(T->CurrArray + startx*grid_size, grid_size, MPI_DOUBLE, world_rank-1, startx, MPI_COMM_WORLD);
        }
        if (world_rank < world_size - 1)  {  // not last process
            //Send my column "end" to the right process world_rank+1
            MPI_Send(T->CurrArray + endx*grid_size, grid_size, MPI_DOUBLE, world_rank+1, endx, MPI_COMM_WORLD);
            //Receive column "end+1" from the right process world_rank+1, populate local Curr Array
            MPI_Recv(T->CurrArray +((endx+1) * grid_size), grid_size, MPI_DOUBLE, world_rank+1, endx+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }  // odd rank
    // --- synchronization end -----
    T->SwapArrays();

}    // end for local_stepcount
  double time_taken = t1.stop();

  if(world_rank == 0){
       // print process statistics and other results
       printf("%d, %d, %d, %g\n", world_rank, startx, endx, time_taken);

       printf("Temp[%d,%d]=%g\n",0,0,T->temp(0,0));
       for(int i= 1; i <5; i++){
         int index = i*(grid_size/5);
         if(index <= endx && index >=startx){
           printf("Temp[%d,%d]=%g\n",index,index,T->temp(index,index));
         }
       }
       printf("Temp[%d,%d]=%g\n",endx,endx,T->temp(endx,endx));

  }
  else{
      // print process statistics and relevant point temperatures
      printf("%d, %d, %d, %g\n", world_rank, startx, endx, time_taken);
      for(int i= 1; i <5; i++){
        int index = i*(grid_size/5);
        if(index <= endx && index >=startx){
          printf("Temp[%d,%d]=%g\n",index,index,T->temp(index,index));
        }
      }
      printf("Temp[%d,%d]=%g\n",endx,endx,T->temp(endx,endx));
   }


  delete T;
  if(world_rank == 0){
    double global_time = global_timer.stop();
    printf("Time taken (in seconds) : %g\n", global_time);
  }
  MPI_Finalize();
  return 0;
}
