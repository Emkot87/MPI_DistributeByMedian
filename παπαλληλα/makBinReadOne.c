#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


float **alloc_2d_init(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int main(int argc, char* argv[]){
    srand((unsigned int)time(NULL));
    int N, d;
    N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);

    float *pointsForFile = (float*)malloc(N*d*sizeof(float));
    //pointsForFile = alloc_2d_init(N,d);
    for(int i = 0; i < N ; i++){
        for(int j = 0; j < d ; j++){
            float k = (float)rand()/((float)RAND_MAX/100);
            pointsForFile[i*d+j] = k;
            printf("%f\n",k);
        }
    }
    FILE *fptr;
    
    fptr = fopen("program2.bin","wb+");
    
    fwrite(pointsForFile, N*d,sizeof(float),fptr);
    fclose(fptr);
    
    printf("------------------\n");

    float *vals = (float*)malloc(N*d*sizeof(float));
    FILE *sad;
    sad = fopen("program2.bin","rb");

    fread(vals,sizeof(float),N*d,sad);
    fclose(sad);

    for(int i = 0 ; i<N ;i++){
        for(int j = 0 ; j < d ; j++){
            printf("%f\n",vals[i*d+j]);
        }
    }
	
}