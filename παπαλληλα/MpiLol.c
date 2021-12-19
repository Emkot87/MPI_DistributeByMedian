#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define swap(x, y) { float temp = x; x = y; y = temp; }

// Standard partition process of QuickSort().
// It considers the last element as pivot
// and moves all smaller element to left of
// it and greater elements to right
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

float* compute_average_to(float** points, float* pivot,int elems, int d){
	float *meanPerProc = (float*)malloc(elems*sizeof(float));
	for(int i = 0 ; i < elems ; i++){
		float s = 0;
		for(int k = 0 ; k < d ; k++){
			s += pow(points[i][k]-pivot[k],2);
		}
		meanPerProc[i] = sqrt(s);
		}
	return meanPerProc;
}

float **alloc_2d_init(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int RandRange(int Min, int Max)
{
    int diff = Max-Min;
    return (int) (((double)(diff+1)/RAND_MAX) * rand() + Min);
}


void distributeByMean(int my_id, int num_procs, int length, int d,float** vals, int first, int last){
	float meanQuick;
	float* distances;
	printf("caller %d entered distribute by mean with %d first an %d last\n",my_id,first,last);
	
	MPI_Status status;
	
	
	if ( first == last ){
		printf("kati teleiwse--------------first %d---------last %d----------------------- \n",first,last);
		return;
	}
	
	if( my_id == first ) {
		//prepei na allazei o arxhgos kathe fora mexri na treksw me 2 ?
		// pws tha to kanw auto 
		// MPI_send fash
		
		int pivot = RandRange(0 , length);
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// send to everyone 
			MPI_Send(vals[pivot],d,MPI_FLOAT,i,1,MPI_COMM_WORLD);
		}
		//printf("%d asked people to do work with pivot 1 %f %d \n",my_id,vals[pivot][1],pivot);
		distances = (float*)malloc(length*sizeof(float));
		float* meanD = (float*)malloc(length*num_procs*sizeof(float));
		
		for(int i = 0 ; i < length ; i++){
			float s = 0;
			for(int k = 0 ; k < d ; k++){
				s += pow(vals[i][k]-vals[pivot][k],2);
			}
			meanD[i] = sqrt(s);
			distances[i] = meanD[i];
		}
		
		// recv the distances from the other procs
		
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone 
			MPI_Recv(&meanD[i*length],length,MPI_FLOAT,i,1,MPI_COMM_WORLD,&status);
		}
		
		for( int i = 0 ; i < length*num_procs ; i ++){
			printf("%d oi meses times einai %f\n" ,i, meanD[i]);
		}
		//printf(" to mhkos pou edwsa einai %d\n",length*num_procs);
		float q1 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 + 1);
		float q2 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 );
		meanQuick = (q1+q2)/2;
		
		printf("\n\n+++++++++++++ diamesos einai %f apo q1 %f kai q2 %f ++++++++++++++++ \n \n",meanQuick, q1,q2);
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone 
			MPI_Send(&meanQuick,1,MPI_FLOAT,i,1,MPI_COMM_WORLD);
		}
		
    }

    else if(my_id > first && my_id <= last) {
		
		float* pvt;
		pvt = (float*)malloc(d*sizeof(float));
		
		// recieve from leader
		MPI_Recv(pvt,d,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
		
        printf("lambanw douleia %d %f\n", my_id,pvt[0]);
		
		distances = (float*)malloc(length*sizeof(float));
	
		for(int i = 0 ; i < length ; i++){
			float s = 0;
			for(int k = 0 ; k < d ; k++){
				s += pow(vals[i][k]-pvt[k],2);
			}
			distances[i] = sqrt(s);
			printf("i'm %d and for i %d my mean dist is %f\n", my_id, i , distances[i]);
		}
		
		MPI_Send(distances,length,MPI_FLOAT,first,1,MPI_COMM_WORLD);
		MPI_Recv(&meanQuick,1,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
		// printf("eimai o %d kai ematha oti h diamesos einai %f\n", my_id,meanQuick);
		
    }
	
	// mporw na dokimasw na kanw tis prakseis apo olous mazi kai na xrhsimopoihsw to gather gia na ta mazepsw
	
	int smaller,bigger;
	
	for(int i = 0 ; i < length ; i++){
		if( distances[i] < meanQuick ){
			smaller++;
		}
	}
	bigger = length - smaller;
	
	printf("eimai o %d kai exw %d mikrotera kai %d megalutera\n",my_id,smaller,bigger);
	
	
	
	//printf("eimai to %d kai kanw mpolikh douleia \n",my_id);
	
	
	// prepei edw na perimenoun kathe fora oses sumetexoun 
	// den thumamai pws to kanw auto 
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*
	if(my_id >= first && my_id <= (last+first)/2 )
		distributeByMean(my_id, num_procs, length/2, d, vals, first, (first + last)/2);
	if( my_id >= (first+last)/2 + 1 && my_id <= last)
		distributeByMean(my_id, num_procs, length/2, d, vals, (first+last)/2 + 1 , last);
	*/
	// sto telos stelnw ston prohgoumeno arxhgo to diko mou tmhma ?
	// kai o arxikos arxhgos tha kserei oloklhro ton pinaka kai tha ton epistrefei h kati , h tha einai o dikos tou allagmenos ?
	
}


int main(int argc,char* argv[]){
	srand((unsigned int)time(NULL));
    
	int ierr, my_id, root_process, num_procs; 
    MPI_Status status;

    ierr = MPI_Init(&argc,&argv);
	
	int N, d;
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
	if(argc < 3){
		printf("bale 2 orismata re lulu, arithmo stoixeiwn kai diastash");
		return 0;
	}
	
	//FILE *fptr;
	
	N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);
	N = N/num_procs ;
    

	// prepei na grapsw ena arxeio kai na to diabazw edw mesa o kathenas ta dika toy 
	// mporw na kanw ena allo c gia na diabazw kai auto aplws na ta diabazei kai na kanei ena akoma arxeio 
	
	// prepei ola na paroun ta arxeia tous apo to file pou eftiakse o 0
	
	float* vals = (float*)malloc(sizeof(float)*N*d);
	for( int i = 0 ; i < N ; i++){
		for( int j = 0 ; j < d ; j++){
			vals[i+j]= i+j;
		}	
	}
	
	MPI_File file;
	int access_all = MPI_MODE_RDONLY;
	MPI_File_open(MPI_COMM_WORLD,"program.bin",access_all,MPI_INFO_NULL,&file);
	MPI_File_seek(file,my_id*N*d*sizeof(float),SEEK_SET);
	int k = N*d;
	MPI_File_read(file,vals,N*d,MPI_FLOAT,&status);
	MPI_File_close(&file);
			
	
    // αρχικα to kathena ena diabazei apo kapoio arxeio iso arithmo apo stoixeia 
	// tha kanw kati gia auto ? arxika oti to kathena kanei ena pinaka me osa stoixeia prepei kai auto
	
	//distributeByMean(my_id , num_procs, N, d, vals, 0, num_procs-1);
	
    
    printf("pe gamwto --------------------------------------------------------------------------%d\n",my_id);
		
	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	
	for( int i = 0 ; i < N ; i++){
		for( int j = 0 ; j < d ; j++){
			printf("%f , my_id %d \n" , vals[i+j],my_id);
		}	
	}


    ierr = MPI_Finalize();
	
	
}