/*
 *  sum_openmp.cpp - Demonstrates parallelism via random fill and sum routines
 *                   This program uses OpenMP.
 */

/*---------------------------------------------------------
 *  Parallel Summation 
 *
 *  1. Each thread generates numints random integers (in parallel OpenMP region)
 *  2. Each thread sums his numints random integers (in parallel OpenMP region)
 *  3  One thread sums the partial results.
 *
 *  NOTE: steps 2-3 are repeated as many times as requested (numiterations)
 *---------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

/*==============================================================
 * print_elapsed (prints timing statistics)
 *==============================================================*/
void print_elapsed(const char* desc, struct timeval* start, struct timeval* end, int niters) {

  struct timeval elapsed;
  /* calculate elapsed time */
  if(start->tv_usec > end->tv_usec) {

    end->tv_usec += 1000000;
    end->tv_sec--;
  }
  elapsed.tv_usec = end->tv_usec - start->tv_usec;
  elapsed.tv_sec  = end->tv_sec  - start->tv_sec;

  printf("\n %s total elapsed time = %ld (usec)",
    desc, (elapsed.tv_sec*1000000 + elapsed.tv_usec) / niters);
}

// Functor to sum the numbers
template<typename T>
struct sum_functor {

  // Constructor
  sum_functor() : m_sum(0) {
  }

  void operator() (int& num) {
    m_sum +=num;
  }

  T get_sum() const {
    return m_sum;
  }

  protected:

  T m_sum;
};

int isPowerOfTwo(int x)
{
  while (((x%2) == 0) && x > 1)
    x = x/2;
  if (x == 1)
    return 1;
  else
    return 0;
}

int closestPowerOfTwo(int x)
{
  int n = 1;
  while ((n*2) < x)
  {
    n = n*2;
  }
  return n;
}
/*==============================================================
 *  Main Program (Parallel Summation)
 *==============================================================*/
int main(int argc, char *argv[]) {

  long numints = 0;
  int numiterations = 0;

  char printIO;

  cout<<"Would you like the data to be displayed? (y/n)";
  cin>>printIO;

  vector<int> data;
  vector<long> partial_sums;

  long total_sum = 0;

  struct timeval start, end;   /* gettimeofday stuff */
  struct timezone tzp;

  if( argc < 3) {
    printf("Usage: %s [numints] [numiterations]\n\n", argv[0]);
    exit(1);
  }

  numints       = atol(argv[1]);
  numiterations = atoi(argv[2]);

  printf("\nExecuting %s: nthreads=%d, numints=%d, numiterations=%d\n",
            argv[0], omp_get_max_threads(), numints, numiterations);

  /* Allocate shared memory, enough for each thread to have numints*/
  data.reserve(numints * omp_get_max_threads());

  /* Allocate shared memory for partial_sums */
  partial_sums.resize(omp_get_max_threads());

  /*****************************************************
   * Generate the random ints in parallel              *
   *****************************************************/
  #pragma omp parallel shared(numints,data) 
  {
    int tid;

    /* get the current thread ID in the parallel region */
    tid = omp_get_thread_num();

    srand(tid + time(NULL));    /* Seed rand functions */

    vector<int>::iterator it_cur = data.begin();
    std::advance(it_cur, tid * numints);

    vector<int>::iterator it_end = it_cur;
    std::advance(it_end, numints);

    for(; it_cur != it_end ; ++it_cur) {

      *it_cur = rand()%1000000;
      //*it_cur = 1;
    }
  }

  if (printIO == 'y'){
    for (int i = 0; i < numints; i++)
      cout<<data[i]<<endl;
    cout<<"Data printed"<<endl;
  }

  /*****************************************************
   * Generate the sum of the ints in parallel          *
   * NOTE: Repeated for numiterations                  *
   *****************************************************/
  gettimeofday(&start, &tzp);
  long n = numints;
  int nthr;
  int hi;
  int lo;
  int work;
  int tid;
  int i;
  int j;

  #pragma omp parallel shared(n,nthr,data,partial_sums) private(i,j,tid,work,lo,hi)
  {
   #pragma omp single
   nthr = omp_get_num_threads();
   if (!isPowerOfTwo(nthr))
    nthr = closestPowerOfTwo(nthr);
   tid = omp_get_thread_num();
   work = (n + nthr-1) / nthr;
   lo = work * tid;
   hi = lo + work;
   if (hi > n)
   hi = n;
   for (i = lo+1; i < hi; i++)
    data[i] = data[i] + data[i-1];
   partial_sums[tid] = data[hi-1];
   #pragma omp barrier
   for (j = 1; j < nthr; j = 2*j)
   {
   if (tid >= j)
   partial_sums[tid] = partial_sums[tid] + partial_sums[tid - j];
   #pragma omp barrier
   }
   for (i = lo; i < hi; i++)
   data[i] = data[i] + partial_sums[tid] - data[hi-1];
  } 

  gettimeofday(&end,&tzp);

  if (printIO == 'y'){
    cout<<endl<<"Prefix sums:"<<endl;
    for (int i = 0; i < numints; i++)
      cout<<data[i]<<endl;
  }

  /*****************************************************
   * Output timing results                             *
   *****************************************************/

  print_elapsed("Summation", &start, &end, numiterations);
  printf("\n Total sum = %6ld\n", total_sum);

  return(0);
}
