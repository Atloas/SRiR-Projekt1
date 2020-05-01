#include "mpi.h" 
#include <stdio.h> 

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

int getDataSize(std::string filename);
void readData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector);
void splitData(int myId, int numProcs, int totalDataSize, int* partialDataStarts, int* partialDataEnds);
void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector);
void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataStarts, int* partialDataEnds);
void broadcastData2(int myId, int numProcs, int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataEnds);
void saveData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, int totalDataSize);

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
	int myId = 0, numProcs = 2;
	std::string filename;
	int totalDataSize = 1000, ownDataSize;
	double dt = 60;     //[s]
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

	if (myId == 0)
	{
		totalDataSize = getDataSize(filename);
	}

#ifdef DEBUG
	std::cout << myId << ": totalDataSize broadcast." << std::endl;
#endif
	MPI_Bcast(&totalDataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double* xPosVector = new double[totalDataSize];  //[m]
	double* yPosVector = new double[totalDataSize];  //[m]
	double* zPosVector = new double[totalDataSize];  //[m]
	double* xVelVector = new double[totalDataSize];  //[m/s]
	double* yVelVector = new double[totalDataSize];  //[m/s]
	double* zVelVector = new double[totalDataSize];  //[m/s]
	double* massVector = new double[totalDataSize];  //[kg]

	if (myId == 0)
	{
		readData(filename, xPosVector, yPosVector, zPosVector, xVelVector, yVelVector, zVelVector, massVector);
	}

#ifdef DEBUG
	std::cout << myId << ": Broadcasting initial data." << std::endl;
#endif
	broadcastInitialData(totalDataSize, xPosVector, yPosVector, zPosVector, xVelVector, yVelVector, zVelVector, massVector);
	splitData(myId, numProcs, totalDataSize, partialDataStarts, partialDataEnds);
	ownDataStart = partialDataStarts[myId];
	ownDataEnd = partialDataEnds[myId];
	ownDataSize = ownDataEnd - ownDataStart + 1;

#ifdef DEBUG
	std::cout << myId << ": Data vectors initialized." << std::endl;
#endif

	//Acceleration vectors are local, their indexes are shifted compared to the global vecotrs, ownDataStart -> 0;
	double* xAccelerationVector = new double[ownDataSize];
	double* yAccelerationVector = new double[ownDataSize];
	double* zAccelerationVector = new double[ownDataSize];

	double xPosDiff;
	double yPosDiff;
	double zPosDiff;
	double r2;
	double magnitude;
	double angleH;
	double angleV;

#ifdef DEBUG
	std::cout << myId << ": Starting simulation." << std::endl;
#endif
	int writeCounter = 0;
	for (double t = 0; t < Tmax; t += dt, writeCounter++)
	{
		for (int i = ownDataStart; i < ownDataEnd + 1; i++)
		{
			xAccelerationVector[i - ownDataStart] = 0;
			yAccelerationVector[i - ownDataStart] = 0;
			zAccelerationVector[i - ownDataStart] = 0;

			for (int j = 0; j < totalDataSize; j++)
			{
				if (i == j)
					continue;

				xPosDiff = xPosVector[i] - xPosVector[j];
				yPosDiff = yPosVector[i] - yPosVector[j];
				zPosDiff = zPosVector[i] - zPosVector[j];
				r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2) + pow(zPosDiff, 2);
				magnitude = G*massVector[j] / r2;
				angleH = atan2(yPosDiff, sqrt(pow(zPosDiff, 2) + pow(xPosDiff, 2)));
				angleV = atan2(zPosDiff, xPosDiff);
				xAccelerationVector[0] += -magnitude*cos(angleH)*cos(angleV);
				yAccelerationVector[0] += -magnitude*sin(angleH);
				zAccelerationVector[0] += -magnitude*cos(angleH)*sin(angleV);
			}
		}

		for (int i = ownDataStart; i < ownDataEnd + 1; i++)
		{
			xVelVector[i] += xAccelerationVector[i - ownDataStart] * dt;
			yVelVector[i] += yAccelerationVector[i - ownDataStart] * dt;
			zVelVector[i] += zAccelerationVector[i - ownDataStart] * dt;
			xPosVector[i] += xVelVector[i] * dt;
			yPosVector[i] += yVelVector[i] * dt;
			zPosVector[i] += zVelVector[i] * dt;
		}

		if (myId == 0 && writeCounter % 10 == 0) {
			saveData(filename, xPosVector, yPosVector, zPosVector, totalDataSize);
		}
		writeCounter++;

		broadcastData(myId, numProcs, xPosVector, yPosVector, zPosVector, partialDataStarts, partialDataEnds);
		//broadcastData2(myId, numProcs, totalDataSize, xPosVector, yPosVector, zPosVector, partialDataEnds);
	}

	delete[] xPosVector;
	delete[] yPosVector;
	delete[] zPosVector;
	delete[] xVelVector;
	delete[] yVelVector;
	delete[] zVelVector;
	delete[] massVector;
	delete[] partialDataStarts;
	delete[] partialDataEnds;

	delete[] xAccelerationVector;
	delete[] yAccelerationVector;
	delete[] zAccelerationVector;

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
	if (datafile.is_open()) {
		while (getline(datafile, line)) {
			++count;
		}
		datafile.close();
	}

#ifdef DEBUG
	std::cout << "0: Data read: " << count - 1 << std::endl;
#endif

	return (count - 1);
}

void splitData(int myId, int numProcs, int totalDataSize, int* partialDataStarts, int* partialDataEnds)
{
#ifdef DEBUG
	std::cout << myId << ": Splitting data." << std::endl;
#endif

	int baseCount = totalDataSize / numProcs;
	int leftover = totalDataSize%numProcs;

	partialDataStarts[0] = 0;
	for (int i = 1; i < numProcs; i++)
	{
		partialDataStarts[i] = partialDataStarts[i - 1] + baseCount;
		if (leftover > 0)
		{
			partialDataStarts[i] += 1;
			leftover--;
		}
		partialDataEnds[i - 1] = partialDataStarts[i] - 1;
	}
	partialDataEnds[numProcs - 1] = totalDataSize - 1;

#ifdef DEBUG
	std::cout << myId << ": Calculated parts:" << std::endl;
	for (int i = 0; i < numProcs; i++)
	{
		std::cout << myId << ": [" << partialDataStarts[i] << ", " << partialDataEnds[i] << "]" << std::endl;
	}
#endif
}

void readData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector)
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
	if (datafile.is_open()) {
		while (getline(datafile, line)) {
			int datapos = 0;
			if (count != -1) {
				while ((pos = line.find(delimiter)) != std::string::npos) {
					std::string token = line.substr(0, pos);
					switch (datapos) {
					case 1:
						xPosVector[count] = atof(token.c_str());
						break;
					case 2:
						yPosVector[count] = atof(token.c_str());
						break;
					case 3:
						zPosVector[count] = atof(token.c_str());
						break;
					case 4:
						xVelVector[count] = atof(token.c_str());
						break;
					case 5:
						yVelVector[count] = atof(token.c_str());
						break;
					case 6:
						zVelVector[count] = atof(token.c_str());
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

void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector)
{
	MPI_Bcast(xPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(yPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(zPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(xVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(yVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(zVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(massVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataStarts, int* partialDataEnds)
{
#ifdef DEBUG
	std::cout << myId << ": Broadcasting partial data." << std::endl;
#endif

	int partialDataSize;
	for (int i = 0; i < numProcs; i++)
	{
		partialDataSize = partialDataEnds[i] - partialDataStarts[i] + 1;
		MPI_Bcast(xPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
		MPI_Bcast(yPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
		MPI_Bcast(zPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
	}

#ifdef DEBUG
	std::cout << myId << ": Finished broadcasting partial data." << std::endl;
#endif
}

void broadcastData2(int myId, int numProcs, int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataEnds)
{
	#ifdef DEBUG
	std::cout << myId << ": Broadcasting partial data." << std::endl;
	#endif
	if (numProcs == 1)
		return;

	MPI_Status status;
	int partialDataSendSize = partialDataEnds[myId] + 1;
	if (myId == 0)
	{
		MPI_Send(xPosVector, partialDataSendSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(yPosVector, partialDataSendSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(zPosVector, partialDataSendSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		int partialDataRecvSize = partialDataEnds[myId - 1] + 1;
		MPI_Recv(xPosVector, partialDataRecvSize, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(yPosVector, partialDataRecvSize, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(zPosVector, partialDataRecvSize, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD, &status);

		if (myId != numProcs - 1)
		{
			MPI_Send(xPosVector, partialDataSendSize, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD);
			MPI_Send(yPosVector, partialDataSendSize, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD);
			MPI_Send(zPosVector, partialDataSendSize, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD);
		}
	}
	
	MPI_Bcast(xPosVector, totalDataSize, MPI_DOUBLE, numProcs - 1, MPI_COMM_WORLD);
	MPI_Bcast(yPosVector, totalDataSize, MPI_DOUBLE, numProcs - 1, MPI_COMM_WORLD);
	MPI_Bcast(zPosVector, totalDataSize, MPI_DOUBLE, numProcs - 1, MPI_COMM_WORLD);

	#ifdef DEBUG
	std::cout << myId << ": Finished broadcasting partial data." << std::endl;
	#endif
}

void saveData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, int totalDataSize)
{
	FILE * resultfile;
	resultfile = fopen("resultdata.txt", "a");

	for (int i = 0; i < totalDataSize; i++)
	{
		fprintf(resultfile, "%d:x;y;z\n", i);
		fprintf(resultfile, "%f;", xPosVector[i]);
		fprintf(resultfile, "%f;", yPosVector[i]);
		fprintf(resultfile, "%f;", zPosVector[i]);
		fprintf(resultfile, "\n");
	}
	fprintf(resultfile, "\n\n");

	fclose(resultfile);
}