/**
 * Ball Simulation
 * 
 *  - Use to simulate and observe an output of how is
 *  ball's velocity at a time and accerolate at a time.
 * 
 *  - Created by Tanatorn Nateesanprasert (Big)
 *    Computer Engineer years 3 59070501035
 * - Created by Patipon Petchtone (Ice)
 *    Computer Engineer years 3 59070501049
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

    /* MPI & OMP Initialize */
    /* Set MPI variable */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Set number of thread */
    numThread = atoi(argv[3]);
    omp_set_num_threads(numThread);

    /* Process Information */
    int processBallStartNumber[numProc];
    int processBallNumber[numProc];
    int processBallEndNumber;
    int tempProcessBallStartNumber;
    int tempProcessBallNumber;

    double startTime;
    double endTime;

    /* Calculated Data */
    int ballNumber;  /* ballNumber - number of balls (N)*/
    int timeStep;    /* timeStep - number of simulation (K)*/
    double interval; /* interval - resolution of time (d)*/

    /* Ball's Information */
    double *ballMass;           /* Array data of all ball's mass */
    double *ballElectricCharge; /* Array data of all ball's charge */
    double *ballXPosition;      /* Array data of all ball's x position */
    double *ballYPosition;      /* Array data of all ball's y position */
    double *ballXVelocity;      /* Array data of all ball's x velocity */
    double *ballYVelocity;      /* Array data of all ball's y velocity */
    int *ballValid;             /* Array data of all ball's status if ball has x or y less than 100000 */
    double *newBallXPosition;   /* Array data of all ball's x new position */
    double *newBallYPosition;   /* Array data of all ball's y new position */
    double *newBallXVelocity;   /* Array data of all ball's x new velocity */
    double *newBallYVelocity;   /* Array data of all ball's y new velocity */
    int *newBallValid;          /* Array data of all ball's status if ball is still valid in the next step */

    /*Force and Accel Information*/
    double force;
    double forceSize;
    double forceSquare;
    double deltaXVector;
    double deltaYVector;
    double cumulativeXForce = 0;
    double cumulativeYForce = 0;

    /* MPI request */
    MPI_Request request[30];

    /* I J K L */
    int i;
    int j;
    int k;
    int l;
    double xi;
    double xj;
    double yi;
    double yj;

    startTime = MPI_Wtime();

    /* MASTER CODE */
    if (id == 0)
    {
        /* Open to read file*/
        fp = fopen(argv[1], "r");
        /* Get ballNumber and timeStep */
        fscanf(fp, "%d %d\n", &ballNumber, &timeStep);
        /* Get interval */
        fscanf(fp, "%lf\n", &interval);
    }

    /* Broadcast data */
    MPI_Bcast(&ballNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timeStep, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&interval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Allocate array for Nth ball */
    ballMass = (double *)malloc(sizeof(double) * ballNumber);
    ballElectricCharge = (double *)malloc(sizeof(double) * ballNumber);
    ballXPosition = (double *)malloc(sizeof(double) * ballNumber);
    ballYPosition = (double *)malloc(sizeof(double) * ballNumber);
    ballXVelocity = (double *)malloc(sizeof(double) * ballNumber);
    ballYVelocity = (double *)malloc(sizeof(double) * ballNumber);
    ballValid = (int *)malloc(sizeof(int) * ballNumber);

    /* Allocate array for next step Nth ball */
    newBallXPosition = (double *)malloc(sizeof(double) * ballNumber);
    newBallYPosition = (double *)malloc(sizeof(double) * ballNumber);
    newBallXVelocity = (double *)malloc(sizeof(double) * ballNumber);
    newBallYVelocity = (double *)malloc(sizeof(double) * ballNumber);
    newBallValid = (int *)malloc(sizeof(int) * ballNumber);

    /* MASTER read file */
    if (id == 0)
    {
        /* Get ball's information */
        for (i = 0; i < ballNumber; i++)
        {
            fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &ballMass[i], &ballElectricCharge[i], &ballXPosition[i], &ballYPosition[i], &ballXVelocity[i], &ballYVelocity[i]);
            /* 1 = still valid and 0 = not valid anymore */
            ballValid[i] = 1;
        }
        /* close input file */
        fclose(fp);

        /* Open the output file and wait for calculated data to write */
        fp = fopen(argv[2], "w");
        for (l = 0; l < ballNumber; l++)
            fprintf(fp, "%lf %lf\r\n", ballXPosition[l], ballYPosition[l]);
        fprintf(fp, "%s\r\n", "---");
    }

    /* Broadcast ball's information */
    MPI_Bcast(ballMass, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballElectricCharge, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballXPosition, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballYPosition, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballXVelocity, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballYVelocity, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ballValid, ballNumber, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        tempProcessBallNumber = ballNumber / numProc;
        int modProcessBallNumber = ballNumber % numProc;
        for (i = 0; i < numProc; i++)
        {
            if (i == 0)
            {
                processBallNumber[i] = tempProcessBallNumber;
                processBallStartNumber[i] = 0;
            }
            else
            {
                processBallNumber[i] = tempProcessBallNumber;
                processBallStartNumber[i] = processBallStartNumber[i - 1] + processBallNumber[i - 1];
                if (modProcessBallNumber > 0)
                {
                    processBallNumber[i]++;
                    modProcessBallNumber--;
                }
            }
        }
    }

    /* broadcast processBallNumber and processBallStartNumber */
    MPI_Bcast(processBallNumber, numProc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(processBallStartNumber, numProc, MPI_INT, 0, MPI_COMM_WORLD);

    /* Calculate processBallEndNumber */
    tempProcessBallStartNumber = processBallStartNumber[id];
    processBallEndNumber = tempProcessBallStartNumber + processBallNumber[id];
    /* Looping K time (Simulation number retreive from input file */
    for (k = 0; k < timeStep; k++)
    {
        /* Looping all ball in each rank to calculate accels and velocity */
#pragma omp parallel for private(i, j, force, deltaXVector, deltaYVector, forceSize, forceSquare) reduction(+ \
                                                                                                            : cumulativeXForce, cumulativeYForce)
        for (i = tempProcessBallStartNumber; i < processBallEndNumber; i++)
        {
            cumulativeXForce = 0;
            cumulativeYForce = 0;
            xi = ballXPosition[i];
            yi = ballYPosition[i];

            if (ballValid[i])
            {
                /* Looping all ball */
                for (j = 0; j < ballNumber; j++)
                {
                    /* Second ball is not equal to the first one */
                    if (ballValid[j] && j != i)
                    {
                        xj = ballXPosition[j];
                        yj = ballYPosition[j];
                        deltaXVector = xi - xj;
                        deltaYVector = yi - yj;
                        forceSquare = (deltaXVector * deltaXVector) + (deltaYVector * deltaYVector);
                        /* Calculate Electric Force on the ball */
                        /* F = (q1 x q2) / ((xi - xj)^2 + (yi - yj)^2) */
                        force = (ballElectricCharge[i] * ballElectricCharge[j]) / (forceSquare);
                        /* Find force vector in x axis */
                        // deltaXVector = xi - xj;
                        /* Find force vector in y axis */
                        // deltaYVector = yi - yj;
                        /* size of force = sqrt((xj-xi)^2 + (yj-yi)^2)) */
                        forceSize = sqrt(forceSquare);
                        /* Calculate cumulativeXForce vector and cumulativeYForce vector */
                        cumulativeXForce += force * deltaXVector / forceSize;
                        cumulativeYForce += force * deltaYVector / forceSize;
                    }
                }
                /* Summarise and calculate ball's accelorate */
                // accelerateX = cumulativeXForce / ballMass[i];
                // accelerateY = cumulativeYForce / ballMass[i];

                /* Collect next ball's position and velocity */
                newBallXPosition[i] = xi + (interval * ballXVelocity[i]);
                newBallYPosition[i] = yi + (interval * ballYVelocity[i]);
                newBallXVelocity[i] = ballXVelocity[i] + (interval * cumulativeXForce / ballMass[i]);
                newBallYVelocity[i] = ballYVelocity[i] + (interval * cumulativeYForce / ballMass[i]);
                // newBallValid[i] = (abs(newBallXPosition[i]) > 1000000) || (abs(newBallYPosition[i]) > 1000000) ? 0 : 1;
                newBallValid[i] = (newBallXPosition[i] > 1000000 || newBallXPosition[i] < -1000000 || newBallYPosition[i] > 1000000 || newBallYPosition[i] < -1000000) ? 0 : 1;
            }
        }

        /* Every ranks broadcast next round data to others */
        for (i = 0; i < numProc; i++)
        {
            tempProcessBallStartNumber = processBallStartNumber[i];
            tempProcessBallNumber = processBallNumber[i];
            MPI_Bcast(&newBallXPosition[tempProcessBallStartNumber], tempProcessBallNumber, MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(&newBallYPosition[tempProcessBallStartNumber], tempProcessBallNumber, MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(&newBallXVelocity[tempProcessBallStartNumber], tempProcessBallNumber, MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(&newBallYVelocity[tempProcessBallStartNumber], tempProcessBallNumber, MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(&newBallValid[tempProcessBallStartNumber], tempProcessBallNumber, MPI_INT, i, MPI_COMM_WORLD);
        }

        /* MASTER write the output file */
        if (id == 0)
        {
            /* write all balls new XPosition, YPosition to the output file */
            for (l = 0; l < ballNumber; l++)
            {
                // printf("rank 0 newValid: %d\n", newBallValid[l]);
                if (newBallValid[l])
                    fprintf(fp, "%lf %lf\r\n", newBallXPosition[l], newBallYPosition[l]);
                else
                    fprintf(fp, "%s\r\n", "out");
            }
            fprintf(fp, "%s\r\n", "---");
        }

        /* Overwritten a prev round XY position and velocity and ball's valid */

        if (k != timeStep - 1)
        {
#pragma omp parallel sections
            {
#pragma omp section
                memcpy(ballXPosition, newBallXPosition, sizeof(double) * ballNumber);
#pragma omp section
                memcpy(ballYPosition, newBallYPosition, sizeof(double) * ballNumber);
#pragma omp section
                memcpy(ballXVelocity, newBallXVelocity, sizeof(double) * ballNumber);
#pragma omp section
                memcpy(ballYVelocity, newBallYVelocity, sizeof(double) * ballNumber);
#pragma omp section
                memcpy(ballValid, newBallValid, sizeof(int) * ballNumber);
            }
        }
    }
    if (id == 0)
    {
        fclose(fp);
        // endTime = MPI_Wtime();
        // printf("total time used : %lf\n", endTime - startTime);
    }
    endTime = MPI_Wtime();
    printf("rank %d use %lf\n", id, endTime - startTime);
    MPI_Finalize();
    return 0;
}
