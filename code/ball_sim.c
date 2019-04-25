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

    /* For Calculate run time */
    double startTime;
    double endTime;

    /* Calculated Data */
    int ballNumber;   /* ballNumber - number of balls (N)*/
    int simNumber;    /* simNumber - number of simulation (K)*/
    double timeStep;  /* timeStep - resolution of time (d)*/
    double totalTime; /* totalTime - total times (T) PS: K = T / d */

    /* Process Information */
    int processBallStartNumber;
    int processBallNumber;
    int processBallBoundNumber;

    /* Ball's Information */
    double *ballMass;             /* Array data of all ball's mass */
    double *ballElectricCharge;   /* Array data of all ball's charge */
    double *ballXPosition;        /* Array data of all ball's x position */
    double *ballYPosition;        /* Array data of all ball's y position */
    double *ballXVelocity;        /* Array data of all ball's x velocity */
    double *ballYVelocity;        /* Array data of all ball's y velocity */
    int *ballValid                /* Array data of all ball's status if ball has x or y less than 100000 */
        double *newBallXPosition; /* Array data of all ball's x new position */
    double *newBallYPosition;     /* Array data of all ball's y new position */
    double *newBallXVelocity;     /* Array data of all ball's x new velocity */
    double *newBallYVelocity;     /* Array data of all ball's y new velocity */
    int *newBallValid             /* Array data of all ball's status if ball is still valid in the next step */

        /*Force and Accel Information*/
        double force;
    double forceSize;
    double deltaXVector;
    double deltaYVector;
    double cumulativeXForce = 0;
    double cumulativeYForce = 0;
    double accelerateX;
    double accelerateY;

    /* I J K */
    int i;
    int j;
    int k;

    /* MPI & OMP Initialize */
    /* Set MPI variable */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Set number of thread */
    numThread = atoi(argv[3]);
    omp_set_num_threads(numThread);

    /* MASTER CODE */
    if (id == 0)
    {
        /* Open to read file*/
        fp = fopen(argv[1], "r");
        /* Get ballNumber and simNumber */
        fscanf(fp, "%d %d\n", &ballNumber, &simNumber);
        /* Get timeStep */
        fscanf(fp, "%lf\n", &timeStep);
    }

    /* Broadcast data */
    MPI_Bcast(&ballNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&simNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    if (id == 0)
    {
        /* Get ball's information */
        for (i = 0; i < ballNumber; i++)
        {
            fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &ballMass[i], &ballElectricCharge[i], &ballXPosition[i], &ballYPosition[i], &ballXVelocity[i], &ballYVelocity[i]);
            /* 1 = still valid and 0 = not valid anymore */
            ballValid[i] = 1;
        }
    }

    /* Broadcast ball's information */
    MPI_Bcast(&ballMass, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballElectricCharge, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballXPosition, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballYPosition, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballXVelocity, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballYVelocity, ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ballValid, ballNumber, MPI_INT, 0, MPI_COMM_WORLD);

    /* Calculate each process ball number */
    processBallNumber = ballNumber / numProc;
    /* Calculate start index of each rank */
    processBallStartNumber = id * processBallNumber;
    /* Calculate each process extras ball except MASTER rank*/
    if (id != 0)
    {
        /* id - 1 is capacity of extras ball that each rank can keep */
        if (ballNumber % numProc <= id - 1)
        {
            /* Increment start index of each rank by extras ball */
            processBallStartNumber = processBallStartNumber + (ballNumber % numProc);
        }
        else
        {
            /* Increment start index of each rank by their max extras ball capacity */
            processBallStartNumber = processBallStartNumber + (id - 1);
        }
        /* Increment ball number of each rank if there are extras ball */
        if (ballNumber % numProc >= id)
        {
            processBallNumber++;
        }
    }

    /* Last ball in each rank */
    processBallBoundNumber = processBallStartNumber + processBallNumber;

    /* Looping K time (Simulation number retreive from input file */
    for (k = 0; k < simNumber; k++)
    {
        /* Looping all ball in each rank to calculate accels and velocity */
        for (i = processBallStartNumber; i < processBallBoundNumber; i++, cumulativeXForce = 0, cumulativeYForce = 0)
        {
            if (ballValid[i])
            {
                /* Looping all ball */
                for (j = 0; j < ballNumber; j++)
                {
                    /* Second ball is not equal to the first one */
                    if (ballValid[j] && j != i)
                    {
                        /* Calculate Electric Force on the ball */
                        /* F = (q1 x q2) / ((xi - xj)^2 + (yi - yj)^2) */
                        force = (ballElectricCharge[i] * ballElectricCharge[j]) / (((ballXPosition[i] - ballXPosition[j]) * (ballXPosition[i] - ballXPosition[j])) + ((ballYPosition[i] - ballYPosition[j]) * (ballYPosition[i] - ballYPosition[j])));

                        /* Find force vector in x axis */
                        deltaXVector = ballXPosition[j] - ballXPosition[i];
                        /* Find force vector in y axis */
                        deltaYVector = ballYPosition[j] - ballYPosition[i];
                        /* size of force = sqrt((xj-xi)^2 + (yj-yi)^2)) */
                        forceSize = sqrt(deltaXVector * deltaXVector + deltaYVector * deltaYVector);
                        /* Calculate cumulativeXForce vector and cumulativeYForce vector */
                        cumulativeXForce += force * deltaXVector / forceSize;
                        cumulativeYForce += force * deltaYVector / forceSize;
                    }
                }
                /* Summarise and calculate ball's accelorate */
                accelerateX = cumulativeXForce / ballMass[i];
                accelerateY = cumulativeYForce / ballMass[i];

                /* Collect next ball's position and velocity */
                newBallXPosition[i] = ballXPosition[i] + (timeStep * ballXVelocity[i]);
                newBallYPosition[i] = ballYPosition[i] + (timeStep * ballYVelocity[i]);
                newBallXVelocity[i] = ballXVelocity[i] + (timeStep * accelerateX);
                newBallYVelocity[i] = ballYVelocity[i] + (timeStep * accelerateY);
                /* If new X or Y position is more than 100000 then this ball is not valid anymore */
                newBallValid[i] = (newBallXPosition[i] > 1000000 || newBallYPosition[i] > 1000000) : 0 ? 1;
            }
        }

        /* Every ranks isend data to adjacent ranks */
    }
}
