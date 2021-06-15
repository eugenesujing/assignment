#include "core/utils.h"
#include "core/cxxopts.h"
#include "core/get_time.h"
#include <pthread.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

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
	double *CurrArray;
	double *PrevArray;
  void assign(double *A, uint x, uint y, double newvalue) {
    A[x*size+y] = newvalue;
  };
  double read(double *A, uint x, uint y) {
    return A[x*size+y];
  }; 
public:	
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
struct argu{
    uint tid;
    uint start;
    uint end;
    uint size;
    uint steps;
    double time;
    TemperatureArray* T;
    CustomBarrier* barrier;
};

 void* heat_transfer_parallel(void* arg) {
     argu* as = (argu*)arg;

  timer t1;
  t1.start();
  uint stepcount;
  for (stepcount = 1; stepcount <= as->steps; stepcount ++) {
	  for (uint x = as->start; x <= as->end; x++) {
		  for (uint y = 0; y < as->size; y++) {
			  as->T->ComputeNewTemp(x, y);
		  }
	  }
  as->barrier->wait(); 
  //pthread_mutex_lock(&mutex);   	
  if (as->tid == 0) {
    //std::cout<<"I'm herer\n";
    // thread 0 should swap arrays. This is the only thread in the serial version
	as->T->SwapArrays();
	as->T->IncrementStepCount();
  as->barrier->wait();
  //pthread_cond_signal(&cond);
  // std::cout<<as->T->ReadStepCount()<<std::endl;
  // std::cout<<stepcount<<std::endl;
	}
	 else {  
	// // other threads should wait until swap is complete
    //  pthread_cond_wait(&cond,&mutex);
    //  pthread_cond_signal(&cond);
    as->barrier->wait();
     //std::cout<<"I'm out\n";
        }
  //       as->T->ReadStepCount();
  //       std::cout<<"------------------------------\n";
  //       std::cout<<as->T->ReadStepCount()<<std::endl;
  // std::cout<<stepcount<<std::endl;
   //pthread_mutex_unlock(&mutex);   
      
	    
    }  // end of current step
    
  as->time = t1.stop();
  pthread_exit(NULL);
}

void heat_transfer_calculation_serial(uint size, uint number_of_threads, TemperatureArray* T, uint steps) {
  timer serial_timer;
  double time_taken = 0.0;
  std::vector<uint> startx(number_of_threads);
  std::vector<uint> endx(number_of_threads);

  // The following code is used to determine start and end of each thread's share of the grid
  // Also used to determine which points to print out at the end of this function
  uint min_columns_for_each_thread = size /   number_of_threads;
  uint excess_columns = size % number_of_threads;
  uint curr_column = 0;

  for (uint i = 0; i < number_of_threads; i++) {
    startx[i] = curr_column;
    if (excess_columns > 0) {
      endx[i] = curr_column + min_columns_for_each_thread;
      excess_columns--;
      } 
    else {
           endx[i] = curr_column + min_columns_for_each_thread - 1;
      }
    curr_column = endx[i]+1;
  }
  CustomBarrier* barrier = new CustomBarrier(number_of_threads); 
  // end of code to determine start and end of each thread's share of the grid
  pthread_t* threads = new pthread_t[number_of_threads];
  argu* args = new argu[number_of_threads];

  serial_timer.start();
  //*------------------------------------------------------------------------

  for(uint j=0; j<number_of_threads;j++){
    args[j].tid=j;
    args[j].start=startx[j];
    args[j].end=endx[j];
    args[j].size=size;
    args[j].steps=steps;
    args[j].T=T;
    args[j].barrier=barrier;
    if(pthread_create(threads+j,NULL,heat_transfer_parallel,args+j)!=0) throw std::runtime_error("failed to start a new thread!");
  }
  // Print these statistics for each thread 
  std::cout << "thread_id, start_column, end_column, time_taken\n";
  for(int j=0;j<number_of_threads;j++){
      if(pthread_join(*(threads+j),NULL)!=0) throw std::runtime_error("failed to join a thread!");
      std::cout << args[j].tid<<", "<<args[j].start<<", " << args[j].end << ", " << std::setprecision(TIME_PRECISION)
              << args[j].time<< "\n";
  }

  delete[] threads;
  delete[] args;
  delete barrier;

  
  uint step = size/5;
  uint position = 0;
  for (uint x = 0; x < 5; x++) {
      std::cout << "Temp[" << position << "," << position << "]=" << T->temp(position,position) << std::endl;
      position += step;
  } 
  // Print temparature at select boundary points;
  for (uint i = 0; i < number_of_threads; i++) {
      std::cout << "Temp[" << endx[i] << "," << endx[i] << "]=" << T->temp(endx[i],endx[i]) << std::endl;
  }

  //*------------------------------------------------------------------------
  time_taken = serial_timer.stop();

  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Heat_transfer_calculation",
                           "Model heat transfer in a grid using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nThreads", "Number of threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
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
  uint n_threads = cl_options["nThreads"].as<uint>();
   if (n_threads < 1) {
	  std::cout << "Number of threads should be equal or greater than 1. Terminating..." << std::endl;
	  return 1;
  }
  uint grid_size = cl_options["gSize"].as<uint>();
   if (grid_size <=0) {
	std::cout << "Grid_size should be positive. Terminating..." << std::endl;
	return 1;
  }
  if(n_threads>grid_size){
    std::cout<<"Grid size should be greater than number of threads. Terminating...\n";
    return 1;
  }
  double init_temp = cl_options["mTemp"].as<double>();
  double Cx = cl_options["iCX"].as<double>();
  double Cy = cl_options["iCY"].as<double>();
  uint steps = cl_options["tSteps"].as<uint>();
   if (steps < 1) {
	std::cout << "Number of threads should be equal or greater than 1. Terminating..." << std::endl;
	return 1;
  }
  std::cout << "Grid Size : " << grid_size << "x" << grid_size << std::endl;
  std::cout << "Number of threads : " << n_threads << std::endl;
  std::cout << "Cx : " << Cx << std::endl << "Cy : " << Cy << std::endl;
  std::cout << "Temperature in the middle of grid : " << init_temp << std::endl;
  std::cout << "Time Steps : " << steps << std::endl;
  
  std::cout << "Initializing Temperature Array..." << std::endl;
  TemperatureArray *T = new TemperatureArray(grid_size, Cx, Cy, init_temp);
  // std::cout<<"Finishing initializing array...\n";
  if (!T) {
      std::cout << "Cannot Initialize Temperature Array...Terminating" << std::endl;
      return 2;
  }
  heat_transfer_calculation_serial (grid_size, n_threads, T, steps);

  delete T;
  return 0;
}
