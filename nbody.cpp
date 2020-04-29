#include "mpi.h" 
#include <stdio.h> 

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

int getDataSize(std::string filename);
void readData(std::string filename, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector);
void splitData(int myId, int numProcs, int totalDataSize, int* ownDataSize, int* partialDataStarts, int* partialDataEnds);
void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector);
void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, int* partialDataStarts, int* partialDataEnds);

/***********************************************
 * TODO:
 * -Expand dataset
 * -Test if data is transferred correctly both in initial (should be fine) and in looped broadcasts (not sure)
 * -Test if it all works, somehow, dunno how to check results
 * -Test result visualization
 * 
 ***********************************************/

int main(int argc, char *argv[]) 
{ 
    int myId = 0, numProcs = 16;
    std::string filename;
    int totalDataSize = 1000, ownDataSize;
    double dt = 3600;     //[s]
    double Tmax = 2.6e6;  //Miesiac
    double G = 6.674e-11;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

#ifdef DEBUG
	std::cout << myId << ": Initialized process" << std::endl;
#endif

    int* partialDataStarts = new int[numProcs];
    int* partialDataEnds = new int[numProcs];
    int ownDataStart, ownDataEnd;

    if(myId == 0)
    {
        totalDataSize = getDataSize(filename);
    }
#ifdef DEBUG
	std::cout << myId << ": totalDataSize broadcast." << std::endl;
#endif
    MPI_Bcast(&totalDataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double* xPosVector = new double[totalDataSize];  //[m]
    double* yPosVector = new double[totalDataSize];  //[m]
    double* xVelVector = new double[totalDataSize];  //[m/s]
    double* yVelVector = new double[totalDataSize];  //[m/s]
    double* massVector = new double[totalDataSize];  //[kg]

    if(myId == 0)
    {
        readData(filename, xPosVector, yPosVector, xVelVector, yVelVector, massVector);
    }
    broadcastInitialData(totalDataSize, xPosVector, yPosVector, xVelVector, yVelVector, massVector);
    splitData(myId, numProcs, totalDataSize, &ownDataSize, partialDataStarts, partialDataEnds);
    ownDataStart = partialDataStarts[myId];
    ownDataEnd = partialDataEnds[myId];

#ifdef DEBUG
	std::cout << myId << ": Data vectors initialized." << std::endl;
#endif

    //Force vectors are local, their indexes are shifted compared to the global vecotrs, ownDataStart -> 0;
    double* xAccelerationVector = new double[ownDataSize];
    double* yAccelerationVector = new double[ownDataSize];

    double xPosDiff;
    double yPosDiff;
    double r2;
    double magnitude;
    double angle;

#ifdef DEBUG
	std::cout << myId << ": Starting simulation." << std::endl;
#endif

    for(double t = 0; t < Tmax; t += dt)
    {
        break;
        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
			xAccelerationVector[i - ownDataStart] = 0;
			yAccelerationVector[i - ownDataStart] = 0;

            for(int j = 0; j < totalDataSize; j++)
            {
                if(i == j)
                    continue;

				//Something is off here
				//Calculated values are fine, but the simulation falls apart??
                xPosDiff = xPosVector[i] - xPosVector[j];
                yPosDiff = yPosVector[i] - yPosVector[j];
                r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2);
                magnitude = G*massVector[j]/r2;
                angle = atan2(yPosDiff, xPosDiff);
                xAccelerationVector[i - ownDataStart] += -magnitude*cos(angle);
                yAccelerationVector[i - ownDataStart] += -magnitude*sin(angle);
            }
        }

        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
            xVelVector[i] += xAccelerationVector[i - ownDataStart]*dt;
            yVelVector[i] += yAccelerationVector[i - ownDataStart]*dt;
            xPosVector[i] += xVelVector[i]*dt;
            yPosVector[i] += yVelVector[i]*dt;
        }

        broadcastData(myId, numProcs, xPosVector, yPosVector, partialDataStarts, partialDataEnds);

        //debugging
        //printf("Time = %f\nEx = %f, Ey = %f\nMx = %f, My = %f\nAngle = %f\n\n", t, xPosVector[0], yPosVector[0], xPosVector[1], yPosVector[1], angle*180/3.1416);
    }

	std::getchar();

    delete[] xPosVector;
    delete[] yPosVector;
    delete[] xVelVector;
    delete[] yVelVector;
    delete[] massVector;
    delete[] partialDataStarts;
    delete[] partialDataEnds;

    delete[] xAccelerationVector;
    delete[] yAccelerationVector;

    MPI_Finalize(); 
    
    return 0; 
} 

int getDataSize(std::string filename)
{
#ifdef DEBUG
	std::cout << "0: Reading data size." << std::endl;
#endif
    //Count number of bodies in datafile
    int count = 0;
    std::string line;
    std::ifstream datafile("nbodydata.txt");
    if(datafile.is_open()){
        while(getline(datafile,line)){
            ++count;     
        }
        datafile.close();
    }
    return (count-1);
}

void splitData(int myId, int numProcs, int totalDataSize, int* ownDataSize, int* partialDataStarts, int* partialDataEnds)
{
#ifdef DEBUG
	std::cout << myId << ": Splitting data." << std::endl;
#endif

    int baseCount = totalDataSize/numProcs;
    int leftover = totalDataSize%numProcs;

    partialDataStarts[0] = 0;
    for(int i = 1; i < numProcs; i++)
    {
        partialDataStarts[i] = partialDataStarts[i - 1] + baseCount;
        if(leftover > 0)
        {
            partialDataStarts[i] += 1;
            leftover--;
        }
        partialDataEnds[i-1] = partialDataStarts[i] - 1;
    }
    partialDataEnds[numProcs - 1] = totalDataSize - 1;
}

void readData(std::string filename, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector)
{
#ifdef DEBUG
	std::cout << "0: Reading file." << std::endl;
#endif
    //Load from datafile
    std::string line;
    std::string delimiter = ";";
    int count = -1;
    int pos;
    std::ifstream datafile("nbodydata.txt");
    if(datafile.is_open()){
        while(getline(datafile,line)){
            int datapos = 0;
            if(count != -1){
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    std::string token = line.substr(0, pos);
                    switch(datapos){
                        case 1:
                            xPosVector[count] = atof(token.c_str());
                            break;
                        case 2:
                            yPosVector[count] = atof(token.c_str());
                            break;
                        case 3:
                            xVelVector[count] = atof(token.c_str());
                            break;
                        case 4:
                            yVelVector[count] = atof(token.c_str());
                            break;
                    }
                    line.erase(0, pos + delimiter.length());
                    datapos++;
                }
                massVector[count] = atof(line.c_str());
            }
            ++count;
        }
        datafile.close();
    }
}

void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector)
{
#ifdef DEBUG
	std::cout << "0: Broadcasting initial data." << std::endl;
#endif

    MPI_Bcast(xPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(massVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, int* partialDataStarts, int* partialDataEnds)
{
#ifdef DEBUG
	std::cout << myId << ": Broadcasting partial data." << std::endl;
#endif

    int partialDataSize;
    for(int i = 0; i < numProcs; i++)
    {
        partialDataSize = partialDataEnds[i] - partialDataStarts[i];
        //Can this work when passing a pointer like this? Will it pull partialDataSize elements starting at the appropriate spot in the array?
        MPI_Bcast(xPosVector+partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(yPosVector+partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}