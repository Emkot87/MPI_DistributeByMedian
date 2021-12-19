#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc,char** argv){
    int N, d;
    int* points[N];
    
    for(int i = 0; i < N; i++){
        points[i] = (int*)malloc(d * sizeof(int));
    }
    


    int ierr, my_id, root_process, num_procs; 
    MPI_Status status;

    ierr = MPI_Init(&argc,&argv);

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // αρχικα το κυριο προσεςς θα στειλει σε ολα τα αλλα τυχαια προσεεσιζ για να εχει το καθενα ισο μεριδιο
    // how do we sent a table tho ?
    //


    if( my_id == 0 ) {
        printf("pe re poust \n");
    }

    else {
        printf("thats some trippy shit %d \n ", my_id);
    }
    
    printf("pe gamwto %d \n",my_id);

    ierr = MPI_Finalize();

}