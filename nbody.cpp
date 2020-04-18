//#include "mpi.h" 
#include <stdio.h> 
#include <math.h> 
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
    
    //Force vectors are local, their indexes are shifted and ownDataStart == 0;
    double* xForceVector = new double[ownDataSize];
    double* yForceVector = new double[ownDataSize];

    if(myid == 0)
    {
        readData(filename, xPosVector, yPosVector, xVelVector, yVelVector, massVector);
    }
    broadcastInitialData(myid, xPosVector, yPosVector, xVelVector, yVelVector, massVector);

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

                double xPosDiff = xPosVector[i] - xPosVector[j];
                double yPosDiff = yPosVector[i] - yPosVector[j];
                double r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2);
                double magnitude = G*massVector[i]*massVector[j]/r2;
                double angle = atan2(yPosDiff, xPosDiff);
                xForceVector[i - ownDataStart] += -magnitude*cos(angle);
                yForceVector[i - ownDataStart] += -magnitude*sin(angle);
            }
        }

        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
            xVelVector[i] += xForceVector[i-ownDataStart]/massVector[i];
            yVelVector[i] += yForceVector[i-ownDataStart]/massVector[i];
            xPosVector[i] += xVelVector[i]*dt;
            yPosVector[i] += yVelVector[i]*dt;
        }

        broadcastData(myid, xPosVector, yPosVector);

        //debugging
        printf("Ex = %f, Ey = %f\nMx = %f, My = %f\n\n", xPosVector[0], yPosVector[0], xPosVector[1], yPosVector[1]);
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