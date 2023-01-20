
// mpic++ -std=c++11 -o solver_4p solver_4p.cpp
// mpirun --hostfile hostfile -np 4 solver_4p


#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "mpi.h"
#define  MASTER   0

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

int heat_equation_cal(int n, int argc, char *argv[]){

  double sum_diff, diff, old_value;
  int numtasks, taskid, task, chunksize, tag1, tag2, tag3, tag4, tag5, tag6, i, j, row, col, source;

 
  double** input_temp = (double **) malloc((n+2) * sizeof(double*));
  double** output_temp = (double **) malloc((n+2) * sizeof(double*));
 
  double* in = (double *) malloc((n+2) * (n+2) * sizeof(double));
  double* out = (double *) malloc((n+2) * (n+2) * sizeof(double));
  for(i = 0; i < n+2; ++i){
    input_temp[i] = &(in[(n+2) * i]);  
  }   
  
  for(i = 0; i < n+2; ++i){
    output_temp[i] = &(out[(n+2) * i]);  
  }

 
  for (i =0; i < n+2; ++i){
    for (j=0; j < n+2; ++j){
      if(i <= n & j <= n & i > 0 & j > 0){
        if (i ==1 | i == n | j == 1 | j == n){
          input_temp[i][j] = 100;
          output_temp[i][j] = 100;
        }else{
          input_temp[i][j] = -100;
          output_temp[i][j] = -100;
        }
      }else{
        input_temp[i][j] = 0;
        output_temp[i][j] = 0;
      }
    }
  }

  double size = (n+2) * (n+2);
 
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Status status;

  if (numtasks != 4) {
    printf("Quitting. Number of MPI tasks must be 4.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(0);
   }


	double c = 0.1;
	double ds = 1.0 / (n+1);
	double dt = (ds)*(ds) * 1.0 / (4*c);
  chunksize = n * 2 / numtasks;
  tag1 = 1;
  tag2 = 2;
  tag3 = 3;
  tag4 = 4;
  tag5 = 5;
  tag6 = 6;
 
  auto t1 = Clock::now();  

  if (taskid == MASTER)
  {
    for (int processor = 1; processor <=3; processor++)
    {
      if (processor == 1)
      {
        row = 0;
        col = chunksize;
      }else if (processor == 2)
      {
        row = chunksize;
        col = 0;
      }else if (processor == 3)
      {
        row = chunksize;
        col = chunksize;
      }
      MPI_Send(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD);
      MPI_Send(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD);
    }
  }else 
  {
    source = MASTER;
    MPI_Recv(&row, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&col, 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
  }


 
  for (int t = 0; t < 1000; t++) { 
    if (taskid == MASTER){

      for (int processor = 1; processor <=3; processor++)
      {
        MPI_Recv(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD, &status);
        MPI_Send(&input_temp[row][col], chunksize + 2, MPI_DOUBLE, processor, tag3, MPI_COMM_WORLD); 
        for (i = row + 1; i <= row + chunksize ; i++ )
        {
          MPI_Send(&input_temp[i][col], 1, MPI_DOUBLE, processor, tag4, MPI_COMM_WORLD);   
          MPI_Send(&input_temp[i][col+chunksize + 1], 1, MPI_DOUBLE, processor, tag5, MPI_COMM_WORLD);           
        } 
        MPI_Send(&input_temp[row + chunksize + 1][col], chunksize + 2, MPI_DOUBLE, processor, tag6, MPI_COMM_WORLD); 
      }
      row = 0;
      col = 0;
     
      for (i = row + 1; i < row + chunksize+1; ++ i ) {
        for (j = col + 1; j < col + chunksize+1; ++ j) {
          output_temp[i][j] = input_temp[i][j]+c*dt/(ds*ds)*(input_temp[i+1][j]+input_temp[i-1][j]-4*input_temp[i][j]+input_temp[i][j+1]+input_temp[i][j-1]);
        }
      }
      swap(input_temp, output_temp);

      for (int processor = 1; processor <=3; processor++)
      {
        MPI_Recv(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD, &status);

        MPI_Recv(&input_temp[row+1][col+1], chunksize, MPI_DOUBLE, processor, tag3, MPI_COMM_WORLD, &status);
        for (i = row + 2 ; i < row + chunksize; i++ )
        { 
          MPI_Recv(&input_temp[i][col+1], 1, MPI_DOUBLE, processor, tag4, MPI_COMM_WORLD, &status);
          MPI_Recv(&input_temp[i][col+chunksize], 1, MPI_DOUBLE, processor, tag5, MPI_COMM_WORLD, &status);
        } 
        MPI_Recv(&input_temp[row+ chunksize][col+1], chunksize, MPI_DOUBLE, processor, tag6, MPI_COMM_WORLD, &status);
      }

    } 

    if (taskid > MASTER){

      task = MASTER;
      MPI_Send(&row, 1, MPI_INT, task, tag1, MPI_COMM_WORLD);
      MPI_Send(&col, 1, MPI_INT, task, tag2, MPI_COMM_WORLD);
      MPI_Recv(&input_temp[row][col], chunksize + 2, MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
      for (i = row + 1; i <= row + chunksize; i++ )
      {          
        MPI_Recv(&input_temp[i][col], 1, MPI_DOUBLE, source, tag4, MPI_COMM_WORLD, &status);
        MPI_Recv(&input_temp[i][col+chunksize + 1], 1, MPI_DOUBLE, source, tag5, MPI_COMM_WORLD, &status);
      } 
      MPI_Recv(&input_temp[row + chunksize +1][col], chunksize + 2, MPI_DOUBLE, source, tag6, MPI_COMM_WORLD, &status);

     
      for (i = row + 1; i < row + chunksize + 1; ++ i ) {
        for (j = col + 1; j < col + chunksize + 1; ++ j) {
          output_temp[i][j] = input_temp[i][j]+c*dt/(ds*ds)*(input_temp[i+1][j]+input_temp[i-1][j]-4*input_temp[i][j]+input_temp[i][j+1]+input_temp[i][j-1]);
        }
      }
      swap(input_temp, output_temp);

      task = MASTER;
      MPI_Send(&row, 1, MPI_INT, task, tag1, MPI_COMM_WORLD);
      MPI_Send(&col, 1, MPI_INT, task, tag2, MPI_COMM_WORLD); 
      MPI_Send(&input_temp[row+1][col+1], chunksize, MPI_DOUBLE, task, tag3, MPI_COMM_WORLD);
      for (i = row + 2; i < row + chunksize; i++ )
      {
        
        MPI_Send(&input_temp[i][col+1], 1, MPI_DOUBLE, task, tag4, MPI_COMM_WORLD);
        MPI_Send(&input_temp[i][col+chunksize], 1, MPI_DOUBLE, task, tag5, MPI_COMM_WORLD);
      } 
      MPI_Send(&input_temp[row + chunksize ][col+1], chunksize, MPI_DOUBLE, task, tag6, MPI_COMM_WORLD);
    } 

  } 
 
  if (taskid == MASTER)
  {
    for (int processor = 1; processor <=3; processor++)
    {
      MPI_Recv(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD, &status);
      MPI_Recv(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD, &status);
      for (i = row + 1 ; i < row + chunksize + 1; i++ )
      { 
        MPI_Recv(&input_temp[i][col+1], chunksize, MPI_DOUBLE, processor, tag3, MPI_COMM_WORLD, &status);
        MPI_Recv(&output_temp[i][col+1], chunksize, MPI_DOUBLE, processor, tag4, MPI_COMM_WORLD, &status);    
      } 
    }
    swap(input_temp, output_temp);
  }else 
  {
    task = MASTER;
    MPI_Send(&row, 1, MPI_INT, task, tag1, MPI_COMM_WORLD);
    MPI_Send(&col, 1, MPI_INT, task, tag2, MPI_COMM_WORLD);
    for (i = row + 1; i < row + chunksize + 1; i++ )
    {  
      MPI_Send(&input_temp[i][col+1], chunksize, MPI_DOUBLE, task, tag3, MPI_COMM_WORLD);
      MPI_Send(&output_temp[i][col+1], chunksize, MPI_DOUBLE, task, tag4, MPI_COMM_WORLD);
    } 
    
  }


  auto t2 = Clock::now();

 
  if (taskid == MASTER){

      double timer = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000.0;
      cout << "N = " << n << "; Processor = " << numtasks << ";  Runtime: "
        << timer
        << " seconds" << std::endl;
  }

  delete [] in;
  delete [] out ;
  delete [] input_temp;
  delete [] output_temp;

  return 0; 
}

int main(int argc, char *argv[]){ 

  MPI_Init(&argc, &argv);
  heat_equation_cal(100, argc, argv);
  heat_equation_cal(1000, argc, argv);
  heat_equation_cal(5000, argc, argv);  
  MPI_Finalize(); 

	return 0;
}