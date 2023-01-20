
// g++ -std=c++11 -o solver_1p solver_1p.cpp
// ./solver_1p

#include <iostream>
#include <algorithm>
using namespace std;

typedef std::chrono::high_resolution_clock Clock;

int heat_equation_cal(int n, int niters){

	double c = 0.1;
	double ds = 1.0 / (n+1);
	double dt = (ds)*(ds) * 1.0 / (4*c);

  	double** input_temp = new double*[n+2];
  	for(int i = 0; i < n+2; ++i){
  		input_temp[i] = new double[n+2];  
  	}  	
	double** output_temp = new double*[n+2];
	for(int i = 0; i < n+2; ++i){
		output_temp[i] = new double[n+2];
	}
    	
    for (int i =0; i < n+2; ++i){
    	for (int j=0; j < n+2; ++j){
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

  auto t1 = Clock::now(); 

	for ( int iter = 0; iter < niters ; iter++ ) {	
		for ( int i = 1; i < n +1; ++ i ) {
			for ( int j = 1; j < n +1; ++ j) {
				output_temp[i][j] = input_temp[i][j]+c*dt/(ds*ds)*(input_temp[i+1][j]+input_temp[i-1][j]-4*input_temp[i][j]+input_temp[i][j+1]+input_temp[i][j-1]);
			}
		}
		swap(input_temp, output_temp);
	}
  auto t2 = Clock::now();

  double timer = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000.0;
  cout << "N = " << n << "; Processor = 1" << ";  Runtime: "
  << timer
    << " seconds" << std::endl;

  return 0; 
}

int main(){
	heat_equation_cal(100, 1000);
  heat_equation_cal(1000, 1000);
  heat_equation_cal(5000, 1000);
	return 0;
}