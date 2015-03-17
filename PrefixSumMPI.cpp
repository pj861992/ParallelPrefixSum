/*
 *  sum_mpi.cpp - Demonstrates parallelism via random fill and sum routines.
 *                This program uses MPI.
 */

/*---------------------------------------------------------
 *  Parallel Summation 
 *
 *  1. each processor generates numints random integers (in parallel)
 *  2. each processor sums his numints random integers (in parallel)
 *  3  Time for processor-wise sums
 *  3.1  All the processes send their sum to Processor 0
 *  3.2  Processor 0 receives the local sum from all the other processes.
 *  3.3  Processor 0 adds up the processor-wise sums (sequentially)
 *  3.4  Processor 0 sends the result to all other processors
 *
 *  NOTE: steps 2-3 are repeated as many times as requested (numiterations)
 *---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <ctime>
#include <unistd.h>

using namespace std;


/*==============================================================
 * p_generate_random_ints (processor-wise generation of random ints)
 *==============================================================*/
void p_generate_random_ints(vector<long>& memory, int n) {

  int i;

  /* generate & write this processor's random integers */
  for (i = 0; i < n; ++i) {

    //memory.push_back(rand()%1000);
    memory.push_back(rand()%5);
  }
  return;
}

/*==============================================================
 * p_summation (processor-wise summation of ints)
 *==============================================================*/

// Functor to sum the numbers
struct sum_functor {

  // Constructor
  sum_functor() : m_sum(0) {
  }

  void operator()(int& num) {
    m_sum +=num;
  }

  long get_sum() const {
    return m_sum;
  }

  protected:

  long m_sum;
};

long p_summation(vector<long>& memory, int start, int count) {

  //sum_functor result = std::for_each(memory.begin(), memory.end(), sum_functor());
  //return result.get_sum();
  //cout << "start = "<<start<<" end = "<< start + count - 1 << endl;
  for (unsigned int i = start + 1; i < start + count; i++)
  {
    memory[i] = memory[i] + memory[i-1];

  }
  return 1;
}

long add_received(vector<long>& memory, long temp, int start, int end) {

  //sum_functor result = std::for_each(memory.begin(), memory.end(), sum_functor());
  //return result.get_sum();
  for (unsigned int i = start; i < end; i++)
  {
    memory[i] = memory[i] + temp;

  }
  return 1;
}

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

void print_node(vector<long>& memory, int start, int end)
{
  for (int i = start; i<=end; i++)
  {
    cout << memory[i]<<endl;
  }
}

int calculate_slice_size(int my_proc_id, int total_procs, int total_data_size)
{
  int data_slice_size;
  /* calculate size of data slice current process should use */ 
  data_slice_size = total_data_size / total_procs;
  /* size of data slice of the last process  */
  if (my_proc_id == total_procs -1) {
      data_slice_size = data_slice_size + (total_data_size%total_procs);
    }

  //cout << "Count for "<< my_proc_id << " is " << data_slice_size << endl;
  return data_slice_size;

}

/*==============================================================
 *  Main Program (Parallel Summation)
 *==============================================================*/
int main(int argc, char **argv) {

  int nprocs, numints, numiterations; /* command line args */

  char printIO;

  cout<<"Would you like the data to be displayed? (y/n)";
  cin>>printIO;

  int my_id, iteration;

  long sum;             /* sum of each individual processor */
  long total_sum;       /* Total sum  */

  vector<long> mymemory; /* Vector to store processes numbers        */
  long* buffer;         /* Buffer for inter-processor communication */

  struct timeval gen_start, gen_end; /* gettimeofday stuff */
  struct timeval start, end;         /* gettimeofday stuff */
  struct timezone tzp;

  MPI_Status status;              /* Status variable for MPI operations */

  /*---------------------------------------------------------
   * Initializing the MPI environment
   * "nprocs" copies of this program will be initiated by MPI.
   * All the variables will be private, only the program owner could see its own variables
   * If there must be a inter-procesor communication, the communication must be
   * done explicitly using MPI communication functions.
   *---------------------------------------------------------*/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id); /* Getting the ID for this process */
  int rank;
  rank = my_id;
  /*ofstream myfile;
  string str1 = "file";
  string str2 = to_string(my_id);
  string str3 = str1 + str2;
  str3 = str3 + ".txt";
  myfile.open (str3);*/

  /*---------------------------------------------------------
   *  Read Command Line
   *  - check usage and parse args
   *---------------------------------------------------------*/

  if(argc < 3) {

    if(my_id == 0)
      printf("Usage: %s [numints] [numiterations]\n\n", argv[0]);

    MPI_Finalize();
    exit(1);
  }

  numints       = atoi(argv[1]);
  numiterations = atoi(argv[2]);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* Get number of processors */
  int size;
  if (!isPowerOfTwo(nprocs))
    nprocs = closestPowerOfTwo(nprocs);
  size = nprocs;
  if (my_id >= nprocs)
  {
    MPI_Finalize();
    exit(1);
  }

  if(my_id == 0)
    printf("\nExecuting %s: nprocs=%d, numints=%d, numiterations=%d\n",
            argv[0], nprocs, numints, numiterations);

  /*---------------------------------------------------------
   *  Initialization
   *  - allocate memory for work area structures and work area
   *---------------------------------------------------------*/
  mymemory.reserve(numints);
  buffer = (long *) malloc(sizeof(long));

  if(buffer == NULL) {

    printf("Processor %d - unable to malloc()\n", my_id);
    MPI_Finalize();
    exit(1);
  }


  /* get starting time */
  gettimeofday(&gen_start, &tzp);
  srand(my_id + time(NULL));     
  int data_slice_size = calculate_slice_size(my_id, nprocs, numints);              /* Seed rand functions */
  p_generate_random_ints(mymemory, data_slice_size);  /* random parallel fill */
  gettimeofday(&gen_end, &tzp);

  if(my_id == 0) {
    print_elapsed("Input generated", &gen_start, &gen_end, 1);
    cout << endl;
  }

  int startindex = 0;
  int endindex = startindex + data_slice_size - 1;
  struct timespec tim, tim2;
   tim.tv_sec = 1;
   tim.tv_nsec = 5000;
   int xx = 0;
  if (printIO == 'y'){
  while(xx<nprocs)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_id == xx)
    {
      print_node(mymemory, startindex, endindex);
      nanosleep(&tim , &tim2);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    xx = xx + 1;
  }

  if (my_id == 0)
    cout<<"Input data displayed"<<endl;}
  
    /*int data_slice_size = calculate_slice_size(my_id, nprocs, numints); 
    int startindex = ((numints / nprocs) * my_id);
    int endindex = startindex + data_slice_size - 1;
    /*Perform a parallel prefix summation at each proc*/
    //p_summation(mymemory, startindex, data_slice_size);
  //int startindex = ((numints / nprocs) * my_id);
  
  /*Perform a parallel prefix summation at each proc*/
  p_summation(mymemory, startindex, data_slice_size);
  //cout << endl<<"Size of memory allcated is "<< mymemory.size()<<endl;
  int last = mymemory.size();
  //cout << mymemory[last - 1] <<endl;


  MPI_Barrier(MPI_COMM_WORLD); /* Global barrier */
  gettimeofday(&start, &tzp);
  int key = 1;
  numiterations = 1;
  long temp = 0;
  long receivedSum = 0;
  //int localSize = mymemory.size();
  //cout << endl << "Local size = " << localSize;
  int leaveOut = 1;
  /* repeat for numiterations times */
  for (iteration = 0; iteration < numiterations; iteration++) {
    while (leaveOut <= nprocs/2)
    {
        if(my_id < nprocs - leaveOut)
        {
            
            MPI_Send(&mymemory[endindex], 1, MPI_INT, my_id + leaveOut,0,MPI_COMM_WORLD);
        }
        if (my_id >= leaveOut)
        {
            MPI_Recv(&temp, 1, MPI_INT, my_id - leaveOut,0,MPI_COMM_WORLD,&status);
            receivedSum += temp;
            mymemory[endindex] += temp;
        }
          leaveOut = leaveOut * 2;
        MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  add_received(mymemory, receivedSum, startindex, endindex);
  MPI_Barrier(MPI_COMM_WORLD);

  if(my_id == 0) {

    gettimeofday(&end,&tzp);

    print_elapsed("Summation", &start, &end, numiterations);
    cout<<endl;
    //printf("\n Total sum = %6ld\n", total_sum);
    
  }
  xx = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  if (printIO == 'y'){
  while(xx<nprocs)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_id == xx)
    {
      print_node(mymemory, startindex, endindex);
      nanosleep(&tim , &tim2);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    xx = xx + 1;
  }
}
  
  //myfile.close();

  /*---------------------------------------------------------
   *  Cleanup
   *  - free memory
   *---------------------------------------------------------*/

  /* free memory */
  free(buffer);

  MPI_Finalize();

  return 0;
} /* main() */
