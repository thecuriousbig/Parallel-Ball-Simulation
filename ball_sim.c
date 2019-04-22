#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

int main(int argc, char **argv)
{
    /* Declare Variables */

    /* Reading input and Writing Output */
    FILE *fp;

    /* For MPI Initialize */
    int numProc;
    int id;

    /* For OMP Initialize */
    int numThread;

    /* For Calculate run time */
    double startTime;
    double endTime;

    /* Calculated Data */
    int ballNumber;        /* ballNumber - number of balls (N)*/
    double simNumber;      /* simNumber - number of simulation (K)*/
    double timeResolution; /* timeResolution - resolution of time (d)*/
    double totalTime       /* totalTime - total times (T) PS: K = T / d */

        /* */

        /* MPI & OMP Initialize */

        /* Set number of thread */
        numThread = atoi(argv[3]);
    omp_set_num_threads(numThread);

    /* Set MPI variable */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Reading input from file */
}
