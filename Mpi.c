#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

struct timeval startwtime, endwtime;
double seq_time;

#define swap(x, y) { float temp = x; x = y; y = temp; }

int partition(float arr[], int l, int r)
{
    float x = arr[r];
	int	i = l;
    for (int j = l; j <= r - 1; j++) {
        if (arr[j] <= x) {
            swap(arr[i], arr[j]);
            i++;
        }
    }
    swap(arr[i], arr[r]);
    return i;
}
 
// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
float kthSmallest(float arr[], int l, int r, int k)
{
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {
 
        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, l, r);
 
        // If position is same as k
        if (index - l == k - 1)
            return arr[index];
 
        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallest(arr, l, index - 1, k);
 
        // Else recur for right subarray
        return kthSmallest(arr, index + 1, r,
                            k - index + l - 1);
    }
 
    // If k is more than number of
    // elements in array
    return 12;
}


int RandRange(int Min, int Max)
{	srand((unsigned int)time(NULL));
    int diff = Max-Min;
    return (int) (((double)(diff+1)/RAND_MAX) * rand() + Min);
}


void distributeByMedian(int iter,float* pivot, int my_id, int num_procs, int length, int d, float* vals, int first, int last){
	float meanQuick;
	float* distances;
	MPI_Request request;
	MPI_Status status;
	
	// if it enters on with itself inside it returns
	if ( first == last ){
		return;
	}
	// we assume that the first process every time is the leader
	if( my_id == first ) {
		// if its the first iteration we have to pick a pivot
		if(iter == 0){
			int pivotIndex = RandRange(0 ,length-1);
			for(int i = 0 ; i < d ; i++){
				pivot[i] = vals[pivotIndex*d + i];
			}
			
			for(int i = first + 1 ; i < last + 1 ; i++ ){
				// send pivot to everyone 
				MPI_Send(pivot,d,MPI_FLOAT,i,1,MPI_COMM_WORLD);
			}
		}

		distances = (float*)malloc(length*sizeof(float));
		float* meanD = (float*)malloc(length*num_procs*sizeof(float));
		
		for(int i = 0 ; i < length ; i++){
			float s = 0;
			for(int k = 0 ; k < d ; k++){
				s += pow(vals[i*d+k]-pivot[k],2);
			}
			meanD[i] = s;
			distances[i] = s;
		}
		
		// recv the distances from the other procs
		
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone 
			MPI_Recv(&meanD[(i-first)*length],length,MPI_FLOAT,i,1,MPI_COMM_WORLD,&status);
		}
		
		float q1 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 + 1);
		float q2 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 );
		meanQuick = (q1+q2)/2;
		
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// send the median to everyone from everyone 
			MPI_Send(&meanQuick,1,MPI_FLOAT,i,1,MPI_COMM_WORLD);
		}
		free(meanD);
		
    }

    else if(my_id > first && my_id <= last) {
		if(iter == 0){
			// recieve from leader
			MPI_Recv(pivot,d,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
			
		}

		distances = (float*)malloc(length*sizeof(float));
		//we calculate the distances to the pivot from every point
		for(int i = 0 ; i < length ; i++){
			float s = 0;
			for(int k = 0 ; k < d ; k++){
				s += pow(vals[i*d+k]-pivot[k],2);
			}
			distances[i] = s;
		}
		// send them to the leader
		MPI_Send(distances,length,MPI_FLOAT,first,1,MPI_COMM_WORLD);
		
	    	//recieve the median
		MPI_Recv(&meanQuick,1,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
		
    }
	

	
	int smaller=0,bigger=0;
	// find how many smaller and bigger points we have compared to pivot
	for(int i = 0 ; i < length ; i++){
		if( distances[i] < meanQuick ){
			smaller++;
		}
	}
	bigger = length - smaller;
	
	//we can skip everything if we have 0 points to trade

	if( !(((my_id <= (last+first)/2) && (bigger == 0)) || ((my_id >= (first+last)/2 + 1) && (smaller == 0))) ){

	
	// partition the points of each process
	if(my_id <= (last+first)/2){

		int	i = 0;
    	for (int j = 0; j <= length - 1; j++) {
        	if (distances[j] > meanQuick) {
            	float temp = distances[j];
				distances[j] = distances[i];
				distances[i] = temp;
				for(int k = 0 ; k < d ; k++){
					float temppoint = vals[i*d + k];
					vals[i*d + k] = vals[j*d + k];
					vals[j*d + k] = temppoint;
				}
            	i++;
        	}
    	}

	}

	// partition the points
	if(my_id >= (first+last)/2 + 1){

		
		int	i = 0;
    	for (int j = 0; j <= length-1; j++) {
        	if (distances[j] < meanQuick) {
            	float temp = distances[j];
				distances[j] = distances[i];
				distances[i] = temp;
				for(int k = 0 ; k < d ; k++){
					float temppoint = vals[i*d + k];
					vals[i*d + k] = vals[j*d + k];
					vals[j*d + k] = temppoint;
				}
            	i++;
        	}
    	}
		
	}
	
    }

	if(my_id == first){
		// recieve from everyone what they want to trade

		int* toTrade = (int*)malloc(sizeof(int)*(last-first+1));
		toTrade[0]=bigger;

		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone how many points they have to trade
			MPI_Recv(&toTrade[i-first],1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
		}

		// o arxhgos tha exei 2 counters/pointers kai tha stelnei se olous me poious na kanoun antallages, aytoi tha kanoyn mexri na tous teleiwsoun ta stoixeia an einai o idios aytos tha kanei aplws to trade
		// sToTrade shows the start of those to give smaller points than the median and bToTrade shows the start of those to give bigger points than the median

		int bToTrade = first; 
		while(toTrade[bToTrade - first] == 0 && bToTrade < (first+last)/2 + 1){
			bToTrade++;
		}

		int sToTrade = (first+last)/2 + 1;
		while(toTrade[sToTrade - first] == 0 && sToTrade < last + 1){
			sToTrade++;
		}

		int tempNumToTrade = 0;

		int srStart = 0;

		
		// tha antallaksoun ta perissotera pithana
		
		
		while(sToTrade <=last && bToTrade < ((first+last)/2 + 1)){
			
		// tha antallaksoun ta perissotera pithana
			
			if(toTrade[sToTrade - first]<toTrade[bToTrade - first]){
				tempNumToTrade = toTrade[sToTrade - first];
				toTrade[sToTrade - first] = 0;
				toTrade[bToTrade - first] -= tempNumToTrade; 
			}

			else{
				tempNumToTrade = toTrade[bToTrade - first];
				toTrade[bToTrade - first] = 0;
				toTrade[sToTrade - first] -= tempNumToTrade;
			}
		

			//posa tha antallakseis
			MPI_Send(&tempNumToTrade,1,MPI_INT,sToTrade,1,MPI_COMM_WORLD);
			//me poion tha antallakseis
			MPI_Send(&bToTrade,1,MPI_INT,sToTrade,1,MPI_COMM_WORLD);

			//printf("to bToTrade einai %d kai to sToTrade einai %d\n",bToTrade,sToTrade);

			if(bToTrade == my_id){
				// ama o arxhgos prepei na antallaksei den xreiazetai na to steilei ara as steilw prwta sto allo
				MPI_Send(&vals[srStart*d],tempNumToTrade*d,MPI_FLOAT,sToTrade,1,MPI_COMM_WORLD);
				MPI_Recv(&vals[srStart*d],tempNumToTrade*d,MPI_FLOAT,sToTrade,1,MPI_COMM_WORLD,&status);
				srStart += tempNumToTrade;
			}
			else{
				// alliws prepei na stelei poios allos me megalutera prepei na steilei
				MPI_Send(&tempNumToTrade,1,MPI_INT,bToTrade,1,MPI_COMM_WORLD);
				MPI_Send(&sToTrade,1,MPI_INT,bToTrade,1,MPI_COMM_WORLD);
			}
				// if they run out of points to send the counter gets to the next process
			if(toTrade[sToTrade - first] == 0)
				sToTrade++;

			if(toTrade[bToTrade - first] == 0)
				bToTrade++;
		}
	}

	// enas apo tous duo tha xreiastei na apothukeusei kapou proswrina ta shmeia pou tha antikatastathoun
	// they temporarily save the points to an array

	else if(my_id > first && (my_id <= (last+first)/2) ){
		
		MPI_Send(&bigger,1,MPI_INT,first,1,MPI_COMM_WORLD);

		// which points to send or where to put the recieving points
		int srStart = 0;

		while(bigger != 0 ){
			
			int tradeWith, tradeNums;
			// oso exei stoixeia na antallaksei tha perimenei na labei me poion tha antallaksei
			MPI_Recv(&tradeNums,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);
			// kai posa stoixeia tha antallaksei me auton
			MPI_Recv(&tradeWith,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);
			


			MPI_Send(&vals[srStart*d],tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD);

			MPI_Recv(&vals[srStart*d],tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD,&status);

			// update where to put files and how many they have to trade
			srStart = srStart + tradeNums;
			bigger -= tradeNums;

		}

	}

	else {
		//printf("I'm %d kai stelnw ta mikrotera stoixeia mou\n",my_id);

		MPI_Send(&smaller,1,MPI_INT,first,1,MPI_COMM_WORLD);

		// which points to send or where to put the recieving points
		int srStart = 0 ;

		while(smaller != 0){
			
			int tradeWith, tradeNums;
			// oso exei stoixeia na antallaksei tha perimenei na labei me poion tha antallaksei
			MPI_Recv(&tradeNums,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);
			// kai posa stoixeia tha antallaksei me auton
			MPI_Recv(&tradeWith,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);

			float* tempPoints = (float*)malloc(sizeof(float)*d*tradeNums);

			MPI_Recv(tempPoints,tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD,&status);

			MPI_Send(&vals[srStart*d],tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD);

			for(int l = 0 ; l < d*tradeNums ; l++ ){
				vals[srStart*d + l] = tempPoints[l];
			}
			// update where to put files and how many they have to trade
			free(tempPoints);
			srStart = srStart + tradeNums;
			smaller -= tradeNums;
			

		}

	}

	free(distances);
	
	
	if(my_id >= first && my_id <= (last+first)/2 )
		distributeByMedian(++iter, pivot, my_id, num_procs/2, length, d, vals, first, (first + last)/2);
	if( my_id >= (first+last)/2 + 1 && my_id <= last)
		distributeByMedian(++iter, pivot, my_id, num_procs/2, length, d, vals, (first+last)/2 + 1 , last);
	
	
}


int main(int argc,char* argv[]){
	int ierr, my_id, root_process, num_procs; 
    MPI_Status status;

    ierr = MPI_Init(&argc,&argv);
	
	int N, d;
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
	if(argc < 3){
		printf("use 2 arguments, how many points and what dimension");
		return 0;
	}
	
	
	N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);
	N = N/num_procs ;
    

	
	// Every file reads its appropriate points from the file makeRand.c created
	
	float* vals = (float*)malloc(sizeof(float)*N*d);
	
	MPI_File file;
	int access_all = MPI_MODE_RDONLY;
	MPI_File_open(MPI_COMM_WORLD,"program2.bin",access_all,MPI_INFO_NULL,&file);
	MPI_File_seek(file,my_id*N*d*sizeof(float),MPI_SEEK_SET);
	int k = N*d;
	MPI_File_read(file,vals,N*d,MPI_FLOAT,&status);
	MPI_File_close(&file);
			
	
	// Initialization of pivot
	float *pivot = (float*)malloc(d*sizeof(float));
	gettimeofday (&startwtime, NULL);
	distributeByMedian(0,pivot, my_id, num_procs, N, d, vals, 0, num_procs-1);
	
    
    
	//wait for everyone to return	
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_id == 0){
		gettimeofday (&endwtime, NULL);
		seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("\n\n-=-=-=-=-=-=-=-+++total time %f with %d points with %d dimensions and %d num_procs+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time,N*num_procs,d,num_procs);
		
	}

	// calculate all the distances to find out if our results are correct
	float* distances = (float*)malloc(sizeof(float)*N);
	
	for(int i = 0 ; i < N ; i++){
		float s = 0;
		for(int k = 0 ; k < d ; k++){
			s += pow(vals[i*d+k]-pivot[k],2);
		}
		distances[i] = sqrt(s);
		
	}

	float max = -1;
	float min = 2147483647; //int max
	for(int i = 0 ; i < N ; i++){
		if(distances[i]<min){
			min = distances[i];
		}
		if(distances[i]>max){
			max = distances[i];
		}
	}
	
	printf("\n\n id: %d minimum %f and maximum %f\n\n",my_id,min,max);

	int check = 1;
	if(my_id == 0){
		float* minmaxs = (float*)malloc(2*num_procs*sizeof(float));
		minmaxs[0] = min;
		minmaxs[1] = max;
		for(int l = 1 ; l < num_procs ; l++){ // recieve mins and maxs
			MPI_Recv(&minmaxs[2*l],1,MPI_FLOAT,l,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&minmaxs[2*l+1],1,MPI_FLOAT,l,1,MPI_COMM_WORLD,&status);
		}
		for(int k = 1 ; k < 2*num_procs ; k++){
			if(minmaxs[k] < minmaxs[k-1]){
				printf("[][][][][][][][][][][][][][][Something went wrong here %d][][][][[]][]]][][][][][][][][][]\n",k);
				check = 0;
			}
		}
		if(check)
			printf("\n----------Everything is correctly placed------------\n\n\n");
		else
			printf("Something went wrong with trading\n");

	}
	else{ // the other processes send their mins and maxs
		MPI_Send(&min,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
		MPI_Send(&max,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
	

	free(pivot);
	free(vals);
    ierr = MPI_Finalize();
	
}
