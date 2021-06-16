#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include<vector>

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

void pageRankParallel(int max_iters,
double& time,
Graph&  g,
PageRankType& pr_curr,
PageRankType& pr_next,
int tid,
CustomBarrier& barrier,
int n_threads){

  timer t1;
  uintV n = g.n_;
  int interval = n/n_threads;
  int step = interval;
  if(tid==n_threads-1){
    step += n%n_threads;
  }
  t1.start();

  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'v', process all its inNeighbors 'u'
    for (uintV v = tid*interval; v < (tid*interval + step); v++) {
      uintE in_degree = g.vertices_[v].getInDegree();
      for (uintE i = 0; i < in_degree; i++) {
        uintV u = g.vertices_[v].getInNeighbor(i);
        uintE u_out_degree = g.vertices_[u].getOutDegree();
        if (u_out_degree > 0)
            pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);
      }
    }
    barrier.wait();
    for (uintV v = 0; v < n; v++) {
      pr_next[v] = PAGE_RANK(pr_next[v]);

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      pr_next[v] = 0.0;
    }
    barrier.wait();
  }
  time = t1.stop();
}


void pageRankSerial(Graph &g, int max_iters, int n_threads) {
  uintV n = g.n_;

  std::vector<PageRankType> pr_curr(n,INIT_PAGE_RANK);
  std::vector<PageRankType> pr_next(n,0);
  std::vector<double> times(n_threads);

  std::vector<std::thread> threads;
  // Pull based pagerank
  timer t1;
  double time_taken = 0.0;

  CustomBarrier barrier(n_threads);
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();

  for(int i=0; i<n_threads; i++){
    threads.push_back(std::thread(pageRankParallel,max_iters,times[i],g,pr_curr,pr_next,i,barrier,n_threads));
  }
  std::cout << "thread_id, time_taken" << std::endl;
  for(int j=0; j<n_threads; j++){
    threads[i].join();
    std::cout<<j<<", "<<times[i]<<std::endl;
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
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_pull",
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

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using DOUBLE\n";
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
