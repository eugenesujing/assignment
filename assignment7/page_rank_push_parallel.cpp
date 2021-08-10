#include "core/utils.h"
#include "core/cxxopts.h"
#include "core/get_time.h"
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cstdio>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef double PageRankType;
#endif



int main(int argc, char *argv[]) {
  timer t0;
  t0.start();
  MPI_Init(NULL,NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of Threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
               {"strategy", "Strategy of decomposition and mapping",
                 cxxopts::value<uint>()->default_value(
                     DEFAULT_STRATEGY)},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  uint strategy = cl_options["strategy"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  if(n_threads<1){
    std::cout<<"Number of Threads should be positive\n";
    return 1;
  }
  if(max_iterations<1){
    std::cout<<"Number of iterations should be positive\n";
    return 1;
  }
  if(strategy<1 || strategy>3){
    std::cout<<"Strategy should be within the range [1,3]";
    return 1;
  }
  if(world_rank == 0){
    #ifdef USE_INT
      std::cout << "Using INT" << "\n";
    #else
      std::cout << "Using DOUBLE" << "\n";
    #endif
      std::cout << std::fixed;
      std::cout << "World Size : " << world_size << "\n";
      std::cout << "Iterations: " << max_iterations << "\n";
      std::cout <<"rank, num_edges, communication_time\n";
  }
  Graph g;
  g.readGraphFromBinary<int>(input_file_path);
  //allocate the vertices to process
  uintV count = 0;

  uintE total_assigned = 0;
  uintV *start = new uintV[world_size];
  uintV *end = new uintV[world_size];

  for(int i=0; i <world_size; i++){
    start[i] = count;
    if(i==world_size-1){
      end[i] = g.n_;
    }else{

      while(total_assigned<(i+1)*(g.m_/world_size)){
        total_assigned += g.vertices_[count].getOutDegree();
        count++;
      }
      end[i] = count;
    }
  }


  PageRankType *pr_curr = new PageRankType[g.n_];
  PageRankType *pr_next = new PageRankType[g.n_];
  for (uintV i = 0; i < g.n_; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }
  double communication_time = 0.0;
  uintE edgesProcessed = 0;
  timer t1;
  uintV* count = new uintv[world_size];
  for(int i=0; i<world_size; i++){
    count[i] = end[i] - start[i];
  }
  for(int i=0; i<max_iterations; i++){
    for(uintV j = start[world_rank]; j < end[world_rank]; j++ ){
      uintE out_degree = g.vertices_[j].getOutDegree();
      edgesProcessed += out_degree;
      for (uintE k = 0; k < out_degree; k++) {
        uintV v = g.vertices_[j].getOutNeighbor(k);
        pr_next[v] += (pr_curr[j] / (PageRankType) out_degree);
      }
    }
    //---synchronization phase 1 start---
    t1.start();
    if(strategy == 1){
      if(world_rank == 0){
        PageRankType* temp = new PageRankType[g.n_];
        for(int c=1; c<world_size; c++){

          //receive next page value from other processes
          MPI_Recv(temp, g.n_, MPI_DOUBLE, c, c, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          for(int d=0; d<g.n_; d++){
            pr_next[d] +=temp[d];
          }

        }
        delete[] temp;
        //send the aggreated results to other processes
        for(int c=1; c<world_size; c++){
          MPI_Send(pr_next+start[world_rank], end[world_rank]-start[world_rank], MPI_DOUBLE, c, c, MPI_COMM_WORLD);
        }
      }else{
        MPI_Send(pr_next, g.n_, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);
        MPI_Recv(pr_next+start[world_rank], end[world_rank]-start[world_rank], MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }else if(strategy==2){
      MPI_Reduce(pr_next, pr_next, g.n_, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Scatterv(pr_next, *count, *start, MPI_DOUBLE, pr_next+start[world_rank], count[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }else{
      //strategy 3
      for(int r=0; r < world_size; r++){

        MPI_Reduce(pr_next+start[r], pr_next+start[r], count[r], MPI_DOUBLE, MPI_SUM, r, MPI_COMM_WORLD);
      }

    }

    //---synchronization phase 1 ends---
    communication_time = t1.stop();
    for (uintV v = start[world_rank]; v < end[world_rank]; v++) {
      pr_next[v] = PAGE_RANK(pr_next[v]);
      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
    }
    //reset pr_next to 0 for all vertices
    for (uintV v = 0; v < g.n_; v++) {
      pr_next[v] = 0.0;
    }

    PageRankType local_sum = 0;
    for (uintV v = start[world_rank]; v < end[world_rank]; v++) {
      local_sum += pr_curr[v];
    }
    //---synchronization phase 2 start---
    if(strategy == 1){
      if(world_rank == 0){
        PageRankType global_sum = 0;
        PageRankType temp = 0;
        printf("%d, %d, %f", world_rank, edgesProcessed, communication_time);
        for(int c=1; c<world_size; c++){
          //receive local sum value from other processes
          MPI_Recv(temp, 1, MPI_DOUBLE, c, c, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          global_sum += temp;

        }
        printf("Sum of page rank : %f",global_sum);
        printf("Time taken (in seconds) : %f",t0.stop());
      }else{
        MPI_Send(local_sum, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);
        printf(%d, %d, %f, world_rank, edgesProcessed, communication_time);
      }
    }else{
      printf("%d, %d, %f", world_rank, edgesProcessed, communication_time);
      PageRankType global_sum = 0;
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(world_rank == 0){

        printf("Sum of page rank : %f",global_sum);
        printf("Time taken (in seconds) : %f",t0.stop());
      }
    }

    //---synchronization phase 2 end---
  }
  delete[] start;
  delete[] end;
  delete[] pr_curr;
  delete[] pr_next;
  MPI_Finalize();
  return 0;
}
