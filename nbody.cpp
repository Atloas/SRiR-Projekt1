//#include "mpi.h" 
#include <stdio.h> 
#include <cmath> 
#include <string>

int getDataSize(std::string filename);
void readData(std::string filename, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector);
void splitData(int myid, int numprocs, int dataSize, int* ownDataSize, int* ownDataStart, int* ownDataEnd);
void broadcastInitialData(int myid, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector);
void broadcastData(int myid, double* xPosVector, double* yPosVector);

int main(int argc, char *argv[]) 
{ 
    int myid = 0, numprocs = 1;
    std::string filename;
    int totalDataSize, ownDataSize, ownDataStart, ownDataEnd;
    double dt = 1;  //[s]
    double Tmax = 2.6e6; //Miesiac
    double G = 6.674e-11;

    //MPI_Init(&argc, &argv); 
    //MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
    //MPI_Comm_rank(MPI_COMM_WORLD, &myid); 

    if(myid == 0)
    {
        totalDataSize = getDataSize(filename);
    }
    //MPI_Bcast(&totalDataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    splitData(myid, numprocs, totalDataSize, &ownDataSize, &ownDataStart, &ownDataEnd);

    double* xPosVector = new double[totalDataSize];  //[m]
    double* yPosVector = new double[totalDataSize];  //[m]
    double* xVelVector = new double[totalDataSize];  //[m/s]
    double* yVelVector = new double[totalDataSize];  //[m/s]
    double* massVector = new double[totalDataSize];  //[kg]
    
    //Force vectors are local, their indexes are shifted compared to the global vecotrs, ownDataStart -> 0;
    double* xForceVector = new double[ownDataSize];
    double* yForceVector = new double[ownDataSize];

    if(myid == 0)
    {
        readData(filename, xPosVector, yPosVector, xVelVector, yVelVector, massVector);
    }
    broadcastInitialData(myid, xPosVector, yPosVector, xVelVector, yVelVector, massVector);

    double xPosDiff;
    double yPosDiff;
    double r2;
    double magnitude;
    double angle;

    for(double t = 0; t < Tmax; t += dt)
    {
        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
			xForceVector[i - ownDataStart] = 0;
			yForceVector[i - ownDataStart] = 0;

            for(int j = 0; j < totalDataSize; j++)
            {
                if(i == j)
                    continue;

				//Something is off here
                xPosDiff = xPosVector[i] - xPosVector[j];
                yPosDiff = yPosVector[i] - yPosVector[j];
                r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2);
                magnitude = G*massVector[i]*massVector[j]/r2;
                angle = atan2(yPosDiff, xPosDiff);
                xForceVector[i - ownDataStart] += -magnitude*cos(angle);
                yForceVector[i - ownDataStart] += -magnitude*sin(angle);
            }
        }

        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
            xVelVector[i] += xForceVector[i-ownDataStart]/massVector[i]*dt;
            yVelVector[i] += yForceVector[i-ownDataStart]/massVector[i]*dt;
            xPosVector[i] += xVelVector[i]*dt;
            yPosVector[i] += yVelVector[i]*dt;
        }

        broadcastData(myid, xPosVector, yPosVector);

        //debugging
        printf("Time = %f\nEx = %f, Ey = %f\nMx = %f, My = %f\nAngle = %f\n\n", t, xPosVector[0], yPosVector[0], xPosVector[1], yPosVector[1], angle*180/3.1416);
    }

	std::getchar();

    delete[] xPosVector;
    delete[] yPosVector;
    delete[] xVelVector;
    delete[] yVelVector;
    delete[] massVector;

    delete[] xForceVector;
    delete[] yForceVector;

    //MPI_Finalize(); 
    
    return 0; 
} 

int getDataSize(std::string filename)
{
    //TODO: Read from file
    return 2;
}

void splitData(int myid, int numprocs, int dataSize, int* ownDataSize, int* ownDataStart, int* ownDataEnd)
{
    //TODO: Split data according to myid
    *ownDataSize = dataSize;
    *ownDataStart = 0;
    *ownDataEnd = dataSize - 1;
}

void readData(std::string filename, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector)
{
    //TODO: Read from file

    //Earth
    xPosVector[0] = 0;
    yPosVector[0] = 0;
    xVelVector[0] = 0;
    yVelVector[0] = 0;
    massVector[0] = 6e24;

    //Moon
    xPosVector[1] = 3.84e8;
    yPosVector[1] = 0;
    xVelVector[1] = 0;
    yVelVector[1] = 1e3;
	massVector[1] = 7.3e22;
}

void broadcastInitialData(int myid, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector)
{
    //TODO
    //Each process receives the same data
}

void broadcastData(int myid, double* xPosVector, double* yPosVector)
{
    //TODO
    //Each process overwrites a portion of its data with data received, based on senderid.
    //Just position? The rest isn't necessary for calcualtions, but it depends what data we want as output.

    //Loop through ids and receive Bcasts, special case for own id.
}