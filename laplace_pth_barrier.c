#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define ITER_MAX 3000 // number of maximum iterations
#define CONV_THRESHOLD 1.0e-5f // threshold of convergence

// function declaration
pthread_barrier_t barrier;

// size of each side of the grid
int size, num_threads;

// matrix to be solved
double **grid;

// auxiliary matrix
double **new_grid;

// array to store the errors of each thread
double *t_err;

// return the maximum value
double max(double a, double b){
    if(a > b)
        return a;
    return b;
}

// return the absolute value of a number
double absolute(double num){
    if(num < 0)
        return -1.0 * num;
    return num;
}

// allocate memory for the grid
void allocate_memory(){
    grid = (double **) malloc(size * sizeof(double *));
    new_grid = (double **) malloc(size * sizeof(double *));

    for(int i = 0; i < size; i++){
        grid[i] = (double *) malloc(size * sizeof(double));
        new_grid[i] = (double *) malloc(size * sizeof(double));
    }
}

// initialize the grid
void initialize_grid(){
    // seed for random generator
    srand(10);

    int linf = size / 2;
    int lsup = linf + size / 10;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            // inicializa região de calor no centro do grid
            if ( i>=linf && i < lsup && j>=linf && j<lsup)
                grid[i][j] = 100;
            else
               grid[i][j] = 0;
            new_grid[i][j] = 0.0;
        }
    }
}

// save the grid in a file
void save_grid(){

    char file_name[30];
    sprintf(file_name, "laplace_barrier.txt");

    // save the result
    FILE *file;
    file = fopen(file_name, "w");

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            fprintf(file, "%lf ", grid[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}
//function declaration
void *Jacobi_iteration(void *args);

int main(int argc, char *argv[]) {

    if(argc != 3){
        printf("Usage: ./pi_pth N T\n");
        printf("N: The size of each side of the domain (grid)\n");
        printf("T: Number of threads\n");
        exit(-1);
    }

    // variables to measure execution time
    struct timeval time_start;
    struct timeval time_end;
    
    size = atoi(argv[1]);
    num_threads = atoi(argv[2]);

    // create an array of p_threads
    pthread_t threads[num_threads];

    // store each thread ID
    int t_id[num_threads];

    // allocate momory for the vector of results
    t_err = (double *) malloc(num_threads * sizeof(double));

    // allocate memory to the grid (matrix)
    allocate_memory();

    // set grid initial conditions
    initialize_grid();

    double err = 1.0;
    int iter = 0;

    printf("Jacobi relaxation calculation: %d x %d grid\n", size, size);

    // get the start time
    gettimeofday(&time_start, NULL);

    // initializes the barrier
    pthread_barrier_init(&barrier, NULL, num_threads);

    // Jacobi iteration
    // This loop will end if either the maximum change reaches below a set threshold (convergence)
    // or a fixed number of maximum iterations have completed
    while(err > CONV_THRESHOLD && iter <= ITER_MAX){
        err = 0.0;

        // create the threads
        for(int i = 0; i < num_threads; i++){
            t_id[i] = i;
            pthread_create(&threads[i], NULL, Jacobi_iteration, (void *) &t_id[i]);
        }

        for(int i = 0; i < num_threads; i++){
            pthread_join(threads[i], NULL);
            err = max(err, t_err[i]);
        }

        if(iter % 100 == 0)
            printf("Error of %0.10lf at iteration %d\n", err, iter);

        iter++;
    }

    // get the end time
    gettimeofday(&time_end, NULL);

    double exec_time = (double) (time_end.tv_sec - time_start.tv_sec) +
                       (double) (time_end.tv_usec - time_start.tv_usec) / 1000000.0;

    //save the final grid in file
    save_grid();

    printf("\nKernel executed in %lf seconds with %d iterations and error of %0.10lf\n", exec_time, iter, err);
    
    free(t_err);

    return 0;
}

void *Jacobi_iteration(void *args){

    // thread id
    int id = *(int *) args;
    // calculate the chunk size
    int chunk = size / num_threads;
    
    // calcute begin and end step of the thread
    int begin = id * chunk;
    int end = begin + chunk;

    // the last thread might have more work
    // if num_steps % num_threads != 0
    if( id == num_threads-1 )
        end = size-2;

    double err = 0.0;

    for( int i = begin+1; i <= end; i++) {
        for(int j = 1; j < size-1; j++) {
            new_grid[i][j] = 0.25 * (grid[i][j+1] + grid[i][j-1] +
                grid[i-1][j] + grid[i+1][j]);

            err = max(err, absolute(new_grid[i][j] - grid[i][j]));
            //printf("Erro da tred %d: %lf\n", id, err);
        }
    }

    t_err[id] = err;

    // awaits completion of the other threads
    pthread_barrier_wait(&barrier);

    // copie the next values into the working array for the next iteration
    for( int i = begin+1; i <= end; i++) {
            for( int j = 1; j < size-1; j++)
                grid[i][j] = new_grid[i][j];
        }

}