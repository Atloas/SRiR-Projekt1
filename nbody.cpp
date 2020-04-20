//#include "mpi.h" 
#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 
#include <string>
#include <iostream>
#include <fstream>

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
        break;
        for(int i = ownDataStart; i < ownDataEnd + 1; i++)
        {
			xForceVector[i - ownDataStart] = 0;
			yForceVector[i - ownDataStart] = 0;

            for(int j = 0; j < totalDataSize; j++)
            {
                if(i == j)
                    continue;

				//Something is not right here
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

void splitData(int myid, int numprocs, int dataSize, int* ownDataSize, int* ownDataStart, int* ownDataEnd)
{
    //TODO: Split data according to myid
    *ownDataSize = dataSize;
    *ownDataStart = 0;
    *ownDataEnd = dataSize - 1;
}

void readData(std::string filename, double* xPosVector, double* yPosVector, double* xVelVector, double* yVelVector, double* massVector)
{
    //Load from datafile
    std::string line;
    std::string delimiter = "-";
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