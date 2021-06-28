#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <pthread.h>
#include <atomic>

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
#define DEFAULT_STRATEGY "1"

static std::atomic<uintV> nextV(0);

struct arguments{
  int max_iters;
  double time;
  Graph * g;
  PageRankType* pr_curr;
  PageRankType* pr_next;
  CustomBarrier* barrier;
  uintV start;
  uintV end;
  int strategy;
  int granularity;
  uintE edgesProcessed;
  uintV verticesProcessed;
  double time_b1;
  double time_b2;
  double time_v;
  int tid;
};

uintV getNextVertexToBeProcessed(uintV granularity){
  return nextV.fetch_add(granularity);
}

void* pageRankParallel(void* arg){
  arguments* args = (arguments*) arg;
  timer t1;
  uintV n = args->g->n_;
  args->time_v = 0;
  args->time_b1 = 0;
  args->time_b2 = 0;
  timer t_v;
  timer t_b1;
  timer t_b2;
  t1.start();

  //std::cout<<"n = "<<n<<std::endl;
  //std::cout<<"strategy="<<args->strategy<<std::endl;
  //std::cout<<"gran="<<args->granularity<<std::endl;
  //std::cout<<"interval = "<<interval<<std::endl;
  //std::cout<<"begin iterations tid: "<<args->tid<<std::endl;
  int max_iters = args->max_iters;
  if(args->strategy==1 || args->strategy==2){
    for (int iter = 0; iter < max_iters; iter++) {
      for (uintV v = args->start; v < args->end; v++) {
        //std::cout<<"enter inner loop tid: "<<args->tid<<std::endl;
        uintE in_degree = args->g->vertices_[v].getInDegree();
        args->edgesProcessed += in_degree;
        for (uintE i = 0; i < in_degree; i++) {
          //std::cout<<"enter second inner loop tid: "<<args->tid<<std::endl;
          uintV u = args->g->vertices_[v].getInNeighbor(i);
          uintE u_out_degree = args->g->vertices_[u].getOutDegree();
          if (u_out_degree > 0)
              args->pr_next[v] += (args->pr_curr[u] / (PageRankType) u_out_degree);
        }
      }
      t_b1.start();
      args->barrier->wait();
      args->time_b1 += t_b1.stop();
      for (uintV v = args->start; v < args->end; v++) {
        args->pr_next[v] = PAGE_RANK(args->pr_next[v]);
        args->verticesProcessed ++;
        // reset pr_curr for the next iteration
        args->pr_curr[v] = args->pr_next[v];
        args->pr_next[v] = 0.0;
      }
      t_b2.start();
      args->barrier->wait();
      args->time_b1 += t_b2.stop();

    }
  }else{
    for (int iter = 0; iter < max_iters; iter++) {
      // for each vertex 'v', process all its inNeighbors 'u'
      //std::cout<<"enter iterations tid: "<<args->tid<<std::endl;

        while(1){
          t_v.start();
          uintV v = getNextVertexToBeProcessed(args->granularity);
          args->time_v += t_v.stop();

          if(v>=n) break;
          for(uintV j=0; j<args->granularity && v<n; j++){
            uintE in_degree = args->g->vertices_[v].getInDegree();
            args->edgesProcessed += in_degree;
            for(uintE i = 0; i< in_degree; i++){
              uintV u = args->g->vertices_[v].getInNeighbor(i);
              uintE u_out_degree = args->g->vertices_[u].getOutDegree();
              if (u_out_degree > 0)
                  args->pr_next[v] += (args->pr_curr[u] / (PageRankType) u_out_degree);
            }
            v++;

          }

          //std::cout<<"-------------k="<<k<<"  iter="<<iter<<"  tid="<<args->tid<<std::endl;

        }

        t_b1.start();
        args->barrier->wait();
        args->time_b1 += t_b1.stop();
        
        if(nextV!=0){
      	  nextV=0;
      	  args->barrier->wait();
        }else{
        args->barrier->wait();
        }
        //std::cout<<"after barrier tid: "<<args->tid<<std::endl;

          while(1){
            t_v.start();
            uintV v = getNextVertexToBeProcessed(args->granularity);
            args->time_v += t_v.stop();
            if(v>=n) break;
            for(uintV j=0; j<args->granularity && v<n; j++){
              args->verticesProcessed ++;
              args->pr_next[v] = PAGE_RANK(args->pr_next[v]);

              // reset pr_curr for the next iteration
              args->pr_curr[v] = args->pr_next[v];
              args->pr_next[v] = 0.0;
              v++;
            }
          }
          t_b2.start();
          args->barrier->wait();
          args->time_b1 += t_b2.stop();
          if(nextV!=0){
      	  nextV=0;
      	  args->barrier->wait();
        }else{
        args->barrier->wait();
        }
    }

  }

  //std::cout<<"end iterations tid: "<<args->tid<<std::endl;
  args->time = t1.stop();
  pthread_exit(NULL);
}

void pageRankSerial(Graph &g, int max_iters, int n_threads, int strategy, int granularity) {
  uintV n = g.n_;
  uintE m = g.m_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  pthread_t* pthreads = new pthread_t[n_threads];
  arguments* arrayArg = new arguments[n_threads];
  CustomBarrier barrier(n_threads);

  // Pull based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  //std::cout<<"pthread create begins\n";
  uintE total_assigned = 0;
  uintV count =0;

  uintV interval = n/n_threads;
  for(int i=0; i<n_threads;i++){
    arrayArg[i].max_iters = max_iters;
    arrayArg[i].g = &g;
    arrayArg[i].pr_curr=pr_curr;
    arrayArg[i].pr_next=pr_next;
    arrayArg[i].barrier = &barrier;
    arrayArg[i].strategy = strategy;
    arrayArg[i].edgesProcessed = 0;
    arrayArg[i].verticesProcessed = 0;
    arrayArg[i].tid=i;
    if(strategy ==1){
      uintV step = interval;
      if(i==n_threads-1){
        step+= n%n_threads;
      }
      arrayArg[i].start = interval*i;
      arrayArg[i].end = interval*i +step;
    }
    if(strategy==2){
      arrayArg[i].start = count;
      if(i==n_threads-1){
        arrayArg[i].end = n;
      }else{
        while(total_assigned<(i+1)*(m/n_threads)){
          total_assigned += g.vertices_[count].getInDegree();
          count++;
        }
        arrayArg[i].end = count;
      }


    }
        if(strategy==3){
      arrayArg[i].granularity=1;
    }
    if(strategy ==4){
      arrayArg[i].granularity=granularity;
    }
    //std::cout<<"gran in creating = "<<arrayArg[i].granularity<<std::endl;
    if(pthread_create(pthreads+i,NULL,pageRankParallel,arrayArg+i)!=0) throw std::runtime_error("Fail to create a new thread");
  }
  //std::cout<<"pthread create ends\n";
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  for(int j=0; j<n_threads; j++){
    if(pthread_join(*(pthreads+j),NULL)!=0) throw std::runtime_error("Faile to join a thread");

    std::cout << j <<", "<< arrayArg[j].verticesProcessed <<", "<< arrayArg[j].edgesProcessed<<", "<< arrayArg[j].time_b1<<", "<< arrayArg[j].time_b2<<", "<< arrayArg[j].time_v <<", "<< arrayArg[j].time<< std::endl;
    // Print the above statistics for each thread
    // Example output for 2 threads:
    // thread_id, time_taken
    // 0, 0.12
    // 1, 0.12
  }

  time_taken = t1.stop();
  // -------------------------------------------------------------------


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
          {"strategy", "Strategy of decomposition and mapping",
            cxxopts::value<uint>()->default_value(
                DEFAULT_STRATEGY)},
          {"granularity", "Granularity of decomposition and mapping",
            cxxopts::value<uint>()->default_value(
                DEFAULT_STRATEGY)},
     });

       auto cl_options = options.parse(argc, argv);
       uint n_threads = cl_options["nThreads"].as<uint>();
       uint max_iterations = cl_options["nIterations"].as<uint>();
       std::string input_file_path = cl_options["inputFile"].as<std::string>();
       uint strategy = cl_options["strategy"].as<uint>();
       uint granularity = cl_options["granularity"].as<uint>();
       if(n_threads<1){
         std::cout<<"Number of Threads should be positive\n";
         return 1;
       }
       if(max_iterations<1){
         std::cout<<"Number of iterations should be positive\n";
         return 1;
       }
       if(strategy<1 || strategy>4){
         std::cout<<"Strategy should be within the range [1,4]";
         return 1;
       }
       if(strategy==4 && granularity<1){
         std::cout<<"Granularity should be positive\n";
         return 1;
       }

     #ifdef USE_INT
       std::cout << "Using INT" << std::endl;
     #else
       std::cout << "Using DOUBLE" << std::endl;
     #endif
       std::cout << std::fixed;
       std::cout << "Number of Threads : " << n_threads << std::endl;
       std::cout << "Strategy : " << strategy << std::endl;
       std::cout << "Granularity : " << granularity << std::endl;
       std::cout << "Iterations : " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";


  pageRankSerial(g, max_iterations, n_threads, strategy, granularity);

  return 0;
}
