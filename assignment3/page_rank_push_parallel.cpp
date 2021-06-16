#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <pthread.h>

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



struct arguments{
  int max_iters;
  double time;
  Graph * g;
  PageRankType* pr_curr;
  PageRankType* pr_next;
  int tid;
  CustomBarrier* barrier;
  int n_threads;
  pthread_mutex_t* locks;
};
//remember to exam the input arguments
void* pageRankParallel(void* arg){
  arguments* args = (arguments*) arg;
  timer t1;
  uintV n = args->g->n_;
  int interval = n/args->n_threads;
  int step = interval;
  if(args->tid==args->n_threads-1){
    step += n%args->n_threads;
  }
  t1.start();
  int max_iters = args->max_iters;
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = args->tid*interval; u < args->tid*interval + step; u++) {
      uintE out_degree = args->g->vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = args->g->vertices_[u].getOutNeighbor(i);
        pthread_mutex_lock(args->locks+v);
        args->pr_next[v] += (args->pr_curr[u] / (PageRankType) out_degree);
        pthread_mutex_unlock(args->locks+v);
      }
    }
    args->barrier->wait();
    for (uintV v = args->tid*interval; v < args->tid*interval + step; v++) {
      args->pr_next[v] = PAGE_RANK(args->pr_next[v]);

      // reset pr_curr for the next iteration
      args->pr_curr[v] = args->pr_next[v];
      args->pr_next[v] = 0.0;
    }
    args->barrier->wait();
  }
  args->time = t1.stop();
  pthread_exit(NULL);
}

void pageRankSerial(Graph &g, int max_iters, int n_threads) {
  uintV n = g.n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];

  pthread_t* pthreads = new pthread_t[n_threads];
  arguments* arrayArg = new arguments[n_threads];
  pthread_mutex_t* locks = new pthread_mutex_t[n];

  CustomBarrier barrier(n_threads);

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
    pthread_mutex_init(locks+i,NULL);
  }

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();

  for(int i=0; i<n_threads;i++){
    arrayArg[i].max_iters = max_iters;
    arrayArg[i].g = &g;
    arrayArg[i].pr_curr=pr_curr;
    arrayArg[i].pr_next=pr_next;
    arrayArg[i].tid = i;
    arrayArg[i].barrier = &barrier;
    arrayArg[i].n_threads = n_threads;
    arrayArg[i].locks = locks;
    if(pthread_create(pthreads+i,NULL,pageRankParallel,arrayArg+i)!=0) throw std::runtime_error("Fail to create a new thread");
  }
  //std::cout<<"pthread create ends\n";
  std::cout << "thread_id, time_taken" << std::endl;
  for(int j=0; j<n_threads; j++){
    if(pthread_join(*(pthreads+j),NULL)!=0) throw std::runtime_error("Faile to join a thread");

    std::cout << arrayArg[j].tid <<", "<< arrayArg[j].time << std::endl;
    // Print the above statistics for each thread
    // Example output for 2 threads:
    // thread_id, time_taken
    // 0, 0.12
    // 1, 0.12
  }

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
  delete[] pthreads;
  delete[] arrayArg;
  delete[] locks;

}

int main(int argc, char *argv[]) {
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
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  if(n_threads<1){
    std::cout<<"Number of Threads should be non-negative\n";
    return 1;
  }
  if(nIterations<1){
    std::cout<<"Number of iterations should be non-negative\n";
    return 1;
  }
#ifdef USE_INT
  std::cout << "Using INT" << std::endl;
#else
  std::cout << "Using DOUBLE" << std::endl;
#endif
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Number of Iterations: " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankSerial(g, max_iterations, n_threads);

  return 0;
}
