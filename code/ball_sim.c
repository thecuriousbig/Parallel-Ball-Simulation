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

    /* Declare Variables */

    /* Reading input and Writing Output */
    FILE *fp;

    /* For Calculate run time */
    double startTime;
    double endTime;

    /* Calculated Data */
    int ballNumber;  /* ballNumber - number of balls (N)*/
    int simNumber;   /* simNumber - number of simulation (K)*/
    double timeStep; /* timeStep - resolution of time (d)*/

    /* Process Information */
    int processBallStartNumber[numProc];
    int processBallNumber[numProc];
    int processBallBoundNumber[numProc];
    MPI_Request request[100];

    int processBallStartNumberTmp;
    int processBallNumberTmp;
    int processBallStartNumberTmpLoop;
    int processBallNumberTmpLoop;
    int processBallBoundNumberTmp;

    int diffBallProcessNumber;
    int diffBallBoundProcessNumber;

    int remainWork;
    double start, end;

    /* Ball's Information */
    int baseBallNumber;
    double *ballMass;           /* Array data of all ball's mass */
    double *ballElectricCharge; /* Array data of all ball's charge */
    double *ballXPosition;      /* Array data of all ball's x position */
    double *ballYPosition;      /* Array data of all ball's y position */
    double *ballXVelocity;      /* Array data of all ball's x velocity */
    double *ballYVelocity;      /* Array data of all ball's y velocity */
    int *ballValid;             /* Array data of all ball's status if ball has x or y less than 100000 */

    double *newBallXPosition; /* Array data of all ball's x new position */
    double *newBallYPosition; /* Array data of all ball's y new position */
    double *newBallXVelocity; /* Array data of all ball's x new velocity */
    double *newBallYVelocity; /* Array data of all ball's y new velocity */
    int *newBallValid;        /* Array data of all ball's status if ball is still valid in the next step */

    double *tmpCumXForce; /* Array data of all ball's x cumForce */
    double *tmpCumYForce; /* Array data of all ball's y cumForce */
    double tempTmpCumXForce;
    double tempTmpCumYForce;
    /*Force and Accel Information*/
    double force;
    double forceSquare;
    double forceSize;
    double deltaXVector;
    double deltaYVector;
    double cumulativeXForce = 0;
    double cumulativeYForce = 0;
    double accelerateX;
    double accelerateY;

    /* I J K L*/
    int i;
    int j;
    int k;
    int l;
    start = MPI_Wtime();
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
    MPI_Ibcast(&ballMass[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Ibcast(&ballElectricCharge[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Ibcast(&ballXPosition[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Ibcast(&ballYPosition[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Ibcast(&ballXVelocity[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Ibcast(&ballYVelocity[0], ballNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request[5]);
    MPI_Ibcast(&ballValid[0], ballNumber, MPI_INT, 0, MPI_COMM_WORLD, &request[6]);

    if (id == 0)
    {
        fp = fopen(argv[2], "w");
        for (l = 0; l < ballNumber; l++)
        {
            fprintf(fp, "%lf %lf\r\n", ballXPosition[l], ballYPosition[l]);
        }
        fprintf(fp, "%s\r\n", "---");
    }
    //printf("After first ibcast %d \n", id);
    if (id == 1)
    {
        /* Calculate work of each process */
        baseBallNumber = ballNumber / numProc;
        processBallNumber[0] = baseBallNumber / 1;
        remainWork = (ballNumber % numProc) + (baseBallNumber - processBallNumber[0]);
        processBallStartNumber[0] = 0;
        processBallBoundNumber[0] = processBallNumber[0];
        //	printf("I am 0 start %d end %d",processBallStartNumber[0],processBallBoundNumber[0]);
        for (i = 1; i < numProc; i++)
        {
            processBallNumber[i] = baseBallNumber + (remainWork / (numProc - 1));
            if (remainWork % (numProc - 1) >= i)
            {
                processBallNumber[i]++;
            }
            processBallStartNumber[i] = processBallBoundNumber[i - 1];
            processBallBoundNumber[i] = processBallStartNumber[i] + processBallNumber[i];
            //      printf("I am %d start %d end %d",i,processBallStartNumber[i],processBallBoundNumber[i]);
        }
    }

    /* Broadcast number of work load */

    MPI_Ibcast(processBallStartNumber, numProc, MPI_INT, 1, MPI_COMM_WORLD, &request[7]);
    MPI_Ibcast(processBallBoundNumber, numProc, MPI_INT, 1, MPI_COMM_WORLD, &request[8]);
    MPI_Ibcast(processBallNumber, numProc, MPI_INT, 1, MPI_COMM_WORLD, &request[9]);

    tmpCumXForce = (double *)malloc(sizeof(double) * ballNumber);
    tmpCumYForce = (double *)malloc(sizeof(double) * ballNumber);

    MPI_Waitall(10, request, MPI_STATUS_IGNORE);

    /* Looping K time (Simulation number retreive from input file */
    //printf("this is %d pbNumber = %d \n", id, processBallBoundNumber[id]);
    //printf("before go in loop %d \n", id);
    processBallStartNumberTmp = processBallStartNumber[id];
    processBallBoundNumberTmp = processBallBoundNumber[id];
    processBallNumberTmp = processBallNumber[id];
    diffBallProcessNumber = ballNumber - processBallNumberTmp;
    diffBallBoundProcessNumber = ballNumber - processBallBoundNumberTmp;
    for (k = 0; k < simNumber; k++)
    {

        /* Looping all ball in each rank to calculate accels and velocity */

#pragma omp parallel for schedule(guided) private(i, j, force, deltaXVector, deltaYVector, forceSize, forceSquare) reduction(+ \
                                                                                                                             : cumulativeXForce, cumulativeYForce)
        for (i = processBallStartNumberTmp; i < processBallBoundNumberTmp; i++)
        {
            cumulativeXForce = 0;
            cumulativeYForce = 0;

            if (ballValid[i])
            {
                /* Looping all ball */
                for (j = 0; j < ballNumber; j++)
                {
                    /* Second ball is not equal to the first one */
                    if (ballValid[j] && j != i)
                    {
                        if (j >= processBallStartNumberTmp && j < processBallBoundNumberTmp && k != 0)
                        {
                            j = processBallBoundNumberTmp - 1;
                            cumulativeXForce += tmpCumXForce[i];
                            cumulativeYForce += tmpCumYForce[i];
                        }
                        else
                        {
                            /* Find force vector in x axis */
                            deltaXVector = ballXPosition[i] - ballXPosition[j];
                            /* Find force vector in y axis */
                            deltaYVector = ballYPosition[i] - ballYPosition[j];
                            forceSquare = (deltaXVector * deltaXVector) + (deltaYVector * deltaYVector);
                            /* Calculate Electric Force on the ball */
                            /* F = (q1 x q2) / ((xi - xj)^2 + (yi - yj)^2) */
                            force = (ballElectricCharge[i] * ballElectricCharge[j]) / (forceSquare);

                            /* size of force = sqrt((xj-xi)^2 + (yj-yi)^2)) */
                            forceSize = sqrt(forceSquare);
                            /* Calculate cumulativeXForce vector and cumulativeYForce vector */
                            cumulativeXForce += force * deltaXVector / forceSize;
                            cumulativeYForce += force * deltaYVector / forceSize;
                        }
                    }
                }
                /* Summarise and calculate ball's accelorate */
                //accelerateX = cumulativeXForce / ballMass[i];
                //accelerateY = cumulativeYForce / ballMass[i];
                //printf("new ball x position from rank %d is %lf \n", id, newBallXPosition[i]);
                /* Collect next ball's position and velocity */
                newBallXPosition[i] = ballXPosition[i] + (timeStep * ballXVelocity[i]);
                newBallYPosition[i] = ballYPosition[i] + (timeStep * ballYVelocity[i]);
                newBallXVelocity[i] = ballXVelocity[i] + (timeStep * (cumulativeXForce / ballMass[i]));
                newBallYVelocity[i] = ballYVelocity[i] + (timeStep * (cumulativeYForce / ballMass[i]));
                /* If new X or Y position is more than 100000 then this ball is not valid anymore */
                newBallValid[i] = (newBallXPosition[i] > 1000000 || newBallXPosition[i] < -1000000 || newBallYPosition[i] > 1000000 || newBallYPosition[i] < -1000000) ? 0 : 1;
            }
        }

        //printf("End round %d \n", id);

        for (l = 0; l < numProc; l++)
        {
            processBallStartNumberTmpLoop = processBallStartNumber[l];
            processBallNumberTmpLoop = processBallNumber[l];
            MPI_Ibcast(&newBallXPosition[processBallStartNumberTmpLoop], processBallNumberTmpLoop, MPI_DOUBLE, l, MPI_COMM_WORLD, &request[l * 5 + 0]);
            MPI_Ibcast(&newBallYPosition[processBallStartNumberTmpLoop], processBallNumberTmpLoop, MPI_DOUBLE, l, MPI_COMM_WORLD, &request[l * 5 + 1]);
            MPI_Ibcast(&newBallXVelocity[processBallStartNumberTmpLoop], processBallNumberTmpLoop, MPI_DOUBLE, l, MPI_COMM_WORLD, &request[l * 5 + 2]);
            MPI_Ibcast(&newBallYVelocity[processBallStartNumberTmpLoop], processBallNumberTmpLoop, MPI_DOUBLE, l, MPI_COMM_WORLD, &request[l * 5 + 3]);
            MPI_Ibcast(&newBallValid[processBallStartNumberTmpLoop], processBallNumberTmpLoop, MPI_INT, l, MPI_COMM_WORLD, &request[l * 5 + 4]);
        }

        /*Prepare data for next round*/
        if (k != simNumber - 1)
        {

#pragma omp parallel sections
            {
#pragma omp section
                memcpy(&ballXPosition[processBallStartNumberTmp], &newBallXPosition[processBallStartNumberTmp], sizeof(double) * processBallNumberTmp);
#pragma omp section
                memcpy(&ballYPosition[processBallStartNumberTmp], &newBallYPosition[processBallStartNumberTmp], sizeof(double) * processBallNumberTmp);
#pragma omp section
                memcpy(&ballXVelocity[processBallStartNumberTmp], &newBallXVelocity[processBallStartNumberTmp], sizeof(double) * processBallNumberTmp);
#pragma omp section
                memcpy(&ballYVelocity[processBallStartNumberTmp], &newBallYVelocity[processBallStartNumberTmp], sizeof(double) * processBallNumberTmp);
#pragma omp section
                memcpy(&ballValid[processBallStartNumberTmp], &newBallValid[processBallStartNumberTmp], sizeof(int) * processBallNumberTmp);
            }
        }
        //printf("After last ibcast %d \n", id);
/*Pre calculation before next Round*/
#pragma omp parallel for schedule(guided) private(i, j, force, deltaXVector, deltaYVector, forceSize, forceSquare) reduction(+ \
                                                                                                                             : tempTmpCumXForce, tempTmpCumYForce)
        for (i = processBallStartNumberTmp; i < processBallBoundNumberTmp; i++)
        {
            tempTmpCumXForce = 0;
            tempTmpCumYForce = 0;
            if (ballValid[i])
            {
                /* Looping all ball */
                for (j = processBallStartNumberTmp; j < processBallBoundNumberTmp; j++)
                {
                    /* Second ball is not equal to the first one */
                    if (ballValid[j] && j != i)
                    {
                        /* Find force vector in x axis */
                        deltaXVector = ballXPosition[i] - ballXPosition[j];
                        /* Find force vector in y axis */
                        deltaYVector = ballYPosition[i] - ballYPosition[j];
                        forceSquare = (deltaXVector * deltaXVector) + (deltaYVector * deltaYVector);
                        /* Calculate Electric Force on the ball */
                        /* F = (q1 x q2) / ((xi - xj)^2 + (yi - yj)^2) */
                        force = (ballElectricCharge[i] * ballElectricCharge[j]) / (forceSquare);
                        /* size of force = sqrt((xj-xi)^2 + (yj-yi)^2)) */
                        forceSize = sqrt(forceSquare);
                        /* Calculate cumulativeXForce vector and cumulativeYForce vector */
                        tempTmpCumXForce += force * deltaXVector / forceSize;
                        tempTmpCumYForce += force * deltaYVector / forceSize;
                    }
                }
            }
            tmpCumXForce[i] = tempTmpCumXForce;
            tmpCumYForce[i] = tempTmpCumYForce;
        }

        /*Wait for Ibcast before start new round*/
        MPI_Waitall(5 * numProc, request, MPI_STATUS_IGNORE);

        //printf("After last wait %d \n", id);
        /* MASTER write the output file */
        if (id == 0)
        {
            /* write all balls new XPosition, YPosition to the output file */
            for (l = 0; l < ballNumber; l++)
            {
                // printf("rank 0 newValid: %d\n", newBallValid[l]);
                if (newBallValid[l])
                {
                    // printf("---");
                    // printf("print %lf %lf\r\n", newBallXPosition[l], newBallYPosition[l]);
                    fprintf(fp, "%lf %lf\r\n", newBallXPosition[l], newBallYPosition[l]);
                }
                else
                    fprintf(fp, "%s\r\n", "out");
            }
            fprintf(fp, "%s\r\n", "---");
        }

        /*Copy Remain New positon to old position*/

        if (k != simNumber - 1)
        {
            if (id == 0)
            {
#pragma omp parallel sections
                {
#pragma omp section
                    memcpy(&ballXPosition[processBallBoundNumberTmp], &newBallXPosition[processBallBoundNumberTmp], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballYPosition[processBallBoundNumberTmp], &newBallYPosition[processBallBoundNumberTmp], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballXVelocity[processBallBoundNumberTmp], &newBallXVelocity[processBallBoundNumberTmp], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballYVelocity[processBallBoundNumberTmp], &newBallYVelocity[processBallBoundNumberTmp], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballValid[processBallBoundNumberTmp], &newBallValid[processBallBoundNumberTmp], sizeof(int) * diffBallProcessNumber);
                }
            }
            else if (id == numProc - 1)
            {

#pragma omp parallel sections
                {
#pragma omp section
                    memcpy(&ballXPosition[0], &newBallXPosition[0], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballYPosition[0], &newBallYPosition[0], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballXVelocity[0], &newBallXVelocity[0], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballYVelocity[0], &newBallYVelocity[0], sizeof(double) * diffBallProcessNumber);
#pragma omp section
                    memcpy(&ballValid[0], &newBallValid[0], sizeof(int) * diffBallProcessNumber);
                }
            }
            else
            {
#pragma omp parallel sections
                {
#pragma omp section
                    memcpy(&ballXPosition[0], &newBallXPosition[0], sizeof(double) * processBallStartNumberTmp);
#pragma omp section
                    memcpy(&ballYPosition[0], &newBallYPosition[0], sizeof(double) * processBallStartNumberTmp);
#pragma omp section
                    memcpy(&ballXVelocity[0], &newBallXVelocity[0], sizeof(double) * processBallStartNumberTmp);
#pragma omp section
                    memcpy(&ballYVelocity[0], &newBallYVelocity[0], sizeof(double) * processBallStartNumberTmp);
#pragma omp section
                    memcpy(&ballValid[0], &newBallValid[0], sizeof(int) * processBallStartNumberTmp);

#pragma omp section
                    memcpy(&ballXPosition[processBallBoundNumberTmp], &newBallXPosition[processBallBoundNumberTmp], sizeof(double) * diffBallBoundProcessNumber);
#pragma omp section
                    memcpy(&ballYPosition[processBallBoundNumberTmp], &newBallYPosition[processBallBoundNumberTmp], sizeof(double) * diffBallBoundProcessNumber);
#pragma omp section
                    memcpy(&ballXVelocity[processBallBoundNumberTmp], &newBallXVelocity[processBallBoundNumberTmp], sizeof(double) * diffBallBoundProcessNumber);
#pragma omp section
                    memcpy(&ballYVelocity[processBallBoundNumberTmp], &newBallYVelocity[processBallBoundNumberTmp], sizeof(double) * diffBallBoundProcessNumber);
#pragma omp section
                    memcpy(&ballValid[processBallBoundNumberTmp], &newBallValid[processBallBoundNumberTmp], sizeof(int) * diffBallBoundProcessNumber);
                }
            }
        }
    }
    end = MPI_Wtime();
    printf("rank %dTime is %lf", id, end - start);
    MPI_Finalize();
}
