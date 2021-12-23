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


void distributeByMean(int iter,float* pivot, int my_id, int num_procs, int length, int d, float* vals, int first, int last){
	float meanQuick;
	float* distances;
	//printf("caller %d entered distribute by mean with %d first an %d last on iter %d\n",my_id,first,last,iter);
	MPI_Request request;
	MPI_Status status;
	
	if ( first == last ){
		//printf("kati teleiwse--------------first %d---------last %d----------------------- \n",first,last);
		return;
	}
	
	if( my_id == first ) {
		//prepei na allazei o arxhgos kathe fora mexri na treksw me 2 ?
		// pws tha to kanw auto 
		// MPI_send fash
		if(iter == 0){
			int pivotIndex = RandRange(0 ,length-1);
			for(int i = 0 ; i < d ; i++){
				pivot[i] = vals[pivotIndex*d + i];
				//printf("%f tuxaio index %d to to valpivot tou i\n",vals[0 + i],pivotIndex);
			}
			
			for(int i = first + 1 ; i < last + 1 ; i++ ){
				// send to everyone 
				MPI_Send(pivot,d,MPI_FLOAT,i,1,MPI_COMM_WORLD);
			}
		}

		//printf("%d asked people to do work with pivot 1 %f %d \n",my_id,vals[pivot][1],pivot);
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
		
		for( int i = 0 ; i < length*num_procs ; i ++){
			//printf("%d oi id %d meses times einai %f\n" ,i,my_id, meanD[i]);
		}
		//printf(" to mhkos pou edwsa einai %d\n",length*num_procs);
		float q1 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 + 1);
		float q2 = kthSmallest(meanD,0,length*num_procs-1,length*num_procs/2 );
		meanQuick = (q1+q2)/2;
		
		//printf("\n\n+++++++++++++ diamesos einai %f apo q1 %f kai q2 %f ++++++++++++++++ \n \n",meanQuick, q1,q2);
		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone 
			MPI_Send(&meanQuick,1,MPI_FLOAT,i,1,MPI_COMM_WORLD);
		}
		free(meanD);
		
    }

    else if(my_id > first && my_id <= last) {
		if(iter == 0){
			// recieve from leader
			MPI_Recv(pivot,d,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
			
			//printf("lambanw douleia %d %f\n", my_id,pivot[0]);
		}

		distances = (float*)malloc(length*sizeof(float));

		for(int i = 0 ; i < length ; i++){
			float s = 0;
			for(int k = 0 ; k < d ; k++){
				s += pow(vals[i*d+k]-pivot[k],2);
			}
			distances[i] = s;
			//printf("i'm %d and for i %d my mean dist is %f\n", my_id, i , distances[i]);
		}
		
		MPI_Send(distances,length,MPI_FLOAT,first,1,MPI_COMM_WORLD);

		MPI_Recv(&meanQuick,1,MPI_FLOAT,first,1,MPI_COMM_WORLD,&status);
		// printf("eimai o %d kai ematha oti h diamesos einai %f\n", my_id,meanQuick);
		
    }
	
	// mporw na dokimasw na kanw tis prakseis apo olous mazi kai na xrhsimopoihsw to gather gia na ta mazepsw
	
	int smaller=0,bigger=0;
	//printf("eimai o %d kai to meanQuick einai %f\n",my_id,meanQuick);
	
	for(int i = 0 ; i < length ; i++){
		if( distances[i] < meanQuick ){
			smaller++;
		}
	}
	bigger = length - smaller;
	
	//printf("eimai o %d kai exw %d mikrotera kai %d megalutera\n",my_id,smaller,bigger);

	if( !(((my_id <= (last+first)/2) && (bigger == 0)) || ((my_id >= (first+last)/2 + 1) && (smaller == 0))) ){

	//printf("\n");
	
	if(my_id <= (last+first)/2){
		//printf("eimai o %d kai mphka sto megaluteroi prwtoi \n",my_id);

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


	if(my_id >= (first+last)/2 + 1){

		//printf("eimai o %d kai mphka sto mikroteroi prwtoi \n",my_id);
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

		// NA THUMITHW NA KOITAKSW AUTO POU BRHKA POY XREIAZETAI NA BALW TO i-FIRST STHN ANADROMH THA GINEI MEGAAALH PIPA

		int* toTrade = (int*)malloc(sizeof(int)*(last-first+1));
		toTrade[0]=bigger;
		//printf("%d has %d to trade\n",my_id,toTrade[0]);

		for(int i = first + 1 ; i < last + 1 ; i++ ){
			// recv from everyone 
			//printf("paw na dextw apo ton %d\n",i);
			MPI_Recv(&toTrade[i-first],1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
			//printf("%d has %d to trade\n",i,toTrade[i-first]);
		}
		//printf("bghka apo thn lipsi posa exoun na steiloun\n\n");

		// o arxhgos tha exei 2 counters/pointers kai tha stelnei se olous me poious na kanoun antallages, aytoi tha kanoyn mexri na tous teleiwsoun ta stoixeia an einai o idios aytos tha kanei aplws to trade
		// sToTrade shows the start of those to give small and psajfiajfp

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

		//printf("to bToTrade einai %d kai toTrade[bToTrade] einai %d to sToTrade einai %d kai toTrade[sToTrade] einai %d\n",bToTrade,toTrade[bToTrade],sToTrade,toTrade[sToTrade]);
		// tha antallaksoun ta perissotera pithana
		//printf("toTrade[(first+last)/2-first] einai iso me %d kai toTrade[last-first] = %d \n",toTrade[(first+last)/2-first],toTrade[last-first]);

		// AMA KATI PAEI POLY LATHOS PSAKSE EDW
		
		while(sToTrade <=last && bToTrade < ((first+last)/2 + 1)){
			
		// tha antallaksoun ta perissotera pithana
			
			//printf("uparxoun akoma stoixeia pou prepei na antallagoun\n");
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
			// ara prepei na antallaksoume tempNumTrade stoixeia to btoTrade me to sToTrade
			// printf("Mesa sto id %d WHILE toTrade[(first+last)/2-first] einai iso me %d kai toTrade[last-first] = %d \n",my_id,toTrade[(first+last)/2-first],toTrade[last-first]);
			// printf("Mesa sto id %d WHILE to bToTrade einai %d kai toTrade[bToTrade] einai %d to sToTrade einai %d kai toTrade[sToTrade] einai %d \n",my_id,bToTrade,toTrade[bToTrade],sToTrade,toTrade[sToTrade]);
			
			// printf("tha ginei id %d antallagh %d stoixeiwn\n",my_id,tempNumToTrade);

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

			if(toTrade[sToTrade - first] == 0)
				sToTrade++;

			if(toTrade[bToTrade - first] == 0)
				bToTrade++;
		}
	}

	// enas apo tous duo tha xreiastei na apothukeusei kapou proswrina ta shmeia pou tha antikatastathoun

	else if(my_id > first && (my_id <= (last+first)/2) ){
		//printf("I'm %d kai stelnw ta megalutera stoixeia mou\n",my_id);

		//MPI_Isend(&bigger,1,MPI_INT,first,1,MPI_COMM_WORLD);
		MPI_Send(&bigger,1,MPI_INT,first,1,MPI_COMM_WORLD);

		// which points to send or where to put the recieving points
		int srStart = 0;

		while(bigger != 0 ){
			//printf("%d mphka sthn loupa epeidh exw %d bigger",my_id,bigger);
			int tradeWith, tradeNums;
			// oso exei stoixeia na antallaksei tha perimenei na labei me poion tha antallaksei
			MPI_Recv(&tradeNums,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);
			// kai posa stoixeia tha antallaksei me auton
			MPI_Recv(&tradeWith,1,MPI_INT,first,1,MPI_COMM_WORLD,&status);
			
			//printf("tha antallaksw %d stoixeia me ton %d",tradeNums,tradeWith);


			MPI_Send(&vals[srStart*d],tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD);

			MPI_Recv(&vals[srStart*d],tradeNums*d,MPI_FLOAT,tradeWith,1,MPI_COMM_WORLD,&status);


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
			//printf("%d mphka sthn loupa epeidh exw %d smaller",my_id,smaller);
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

			free(tempPoints);
			srStart = srStart + tradeNums;
			smaller -= tradeNums;
			

		}

	}

	//ypologizw ksana tis apostaseis gia na dw ama eginan ta trades xwris na steilw tis apostaseis

	
	// antallaksa osa eixa mwre
	// prepei edw na perimenoun kathe fora oses sumetexoun 
	// den thumamai pws to kanw auto 
	

	free(distances);
	
	
	if(my_id >= first && my_id <= (last+first)/2 )
		distributeByMean(++iter, pivot, my_id, num_procs/2, length, d, vals, first, (first + last)/2);
	if( my_id >= (first+last)/2 + 1 && my_id <= last)
		distributeByMean(++iter, pivot, my_id, num_procs/2, length, d, vals, (first+last)/2 + 1 , last);
	

	// sto telos stelnw ston prohgoumeno arxhgo to diko mou tmhma ?
	// kai o arxikos arxhgos tha kserei oloklhro ton pinaka kai tha ton epistrefei h kati , h tha einai o dikos tou allagmenos ?
	
}


int main(int argc,char* argv[]){
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
	
	
	N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);
	N = N/num_procs ;
    

	// prepei na grapsw ena arxeio kai na to diabazw edw mesa o kathenas ta dika toy 
	// mporw na kanw ena allo c gia na diabazw kai auto aplws na ta diabazei kai na kanei ena akoma arxeio 
	
	// prepei ola na paroun ta arxeia tous apo to file pou eftiakse o 0
	
	float* vals = (float*)malloc(sizeof(float)*N*d);
	
	MPI_File file;
	int access_all = MPI_MODE_RDONLY;
	MPI_File_open(MPI_COMM_WORLD,"programMnist.bin",access_all,MPI_INFO_NULL,&file);
	MPI_File_seek(file,my_id*N*d*sizeof(float),MPI_SEEK_SET);
	int k = N*d;
	MPI_File_read(file,vals,N*d,MPI_FLOAT,&status);
	MPI_File_close(&file);
			
	
    // αρχικα to kathena ena diabazei apo kapoio arxeio iso arithmo apo stoixeia 
	// tha kanw kati gia auto ? arxika oti to kathena kanei ena pinaka me osa stoixeia prepei kai auto

	float *pivot = (float*)malloc(d*sizeof(float));
	gettimeofday (&startwtime, NULL);
	distributeByMean(0,pivot, my_id, num_procs, N, d, vals, 0, num_procs-1);
	
    
    printf("pe gamwto ------------------------------------------------------------------------%d\n",my_id);
		
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_id == 0){
		gettimeofday (&endwtime, NULL);
		seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("\n\n-=-=-=-=-=-=-=-+++total time %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
		
	}


	float* distances = (float*)malloc(sizeof(float)*N);
	
	for(int i = 0 ; i < N ; i++){
		float s = 0;
		for(int k = 0 ; k < d ; k++){
			s += pow(vals[i*d+k]-pivot[k],2);
		}
		distances[i] = sqrt(s);
		//printf("i'm %d and for i %d my dist from pivot is %f\n", my_id, i , distances[i]);
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

	printf("\n\negw esti %d kai to elaxisto mou einai %f kai to megisto mou %f\n\n",my_id,min,max);

	int check = 1;
	if(my_id == 0){
		float* minmaxs = (float*)malloc(2*num_procs*sizeof(float));
		minmaxs[0] = min;
		minmaxs[1] = max;
		for(int l = 1 ; l < num_procs ; l++){
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
	else{
		MPI_Send(&min,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
		MPI_Send(&max,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
	

	free(pivot);
	free(vals);
    ierr = MPI_Finalize();
	
}