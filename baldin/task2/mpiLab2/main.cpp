#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

void swap(int &a, int &b)
{
	int temp = a;
	a		 = b;
	b		 = temp;
}

int main(int argc, char* argv[])
{
	int ProcNum  = 0;
	int ProcRank = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	if (argc < 4)
	{
		if(ProcRank == 0)
			cout << "NO ENOUGH ARGUMENTS, WILL BE USED SERVICE DATA" << endl;
	}
	
	int* matrix        = nullptr;
	int* sendCount     = nullptr;
	int* resultVectorP = nullptr;  //Parallel   mult
	int* resultVectorC = nullptr;  //Consistent mult
	int n = 5;
	if (argc > 1)
		n = atoi(argv[1]);
	if (n < 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED SIZE = |SIZE|" << endl;
		n = abs(n);
	}
	else if (n == 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED SIZE = 10" << endl;
		n = 10;
	}
	int* multVector   = new int[n];
	int* displs       = new int[ProcNum];
	int* recvCount    = new int[ProcNum];
	int* gatherDispls = new int[ProcNum];
	int border = 0;
	int boolSendCount[2];
	if (ProcRank == 0)
	{
		int Min = 0;
		int Max = 9;
		if (argc > 2)
			Min  = atoi(argv[2]);
		if (argc > 3)
			Max  = atoi(argv[3]);
		if (Max == Min)
		{
			if (ProcRank == 0)
				cout << "ERROR DATA, WILL BE USED MAX = 9, MIN = 0" << endl;
			Max = 9;
			Min = 0;
		}
		else if (Max < Min)
		{
			if (ProcRank == 0)
				cout << "ERROR DATA, MIN AND MAX WILL BE SWAPPED" << endl;
			swap(Min, Max);
		}
		cout << "USING DATA: " << endl << "Size = " << n << endl << "Min = " << Min << endl << "Max = " << Max << endl;

	
		srand(time(NULL));
		matrix = new int[n * n];

		for (int i = 0; i < n; i++)
		{
			multVector[i] = Min + rand() % (Max - Min + 1);
			for (int j = 0; j < n; j++)
				matrix[i * n + j] = Min + rand() % (Max - Min + 1);
		}	

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				cout << matrix[i * n + j] << "   ";
			cout << "    " << multVector[i] << endl;
		}
		cout << endl << "--------------------------------" << endl << endl;
			
		resultVectorC = new int[n];
		for (int i = 0; i < n; i++)
			resultVectorC[i] = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				resultVectorC[i] += matrix[i * n + j] * multVector[j];      
			
		resultVectorP = new int[n];
		for (int i = 0; i < n; i++)
			resultVectorP[i] = 0;
				
		int startNum = n / ProcNum;
		sendCount = new int[ProcNum];
		for (int i = 0; i < ProcNum; i++)
		{
			sendCount[i] = startNum * n;
			displs[i] = 0;
		}
		border = n - startNum * ProcNum;
		for (int i = 0; i < border; i++)
			sendCount[i] += n;

		boolSendCount[0] = sendCount[0];
		boolSendCount[1] = sendCount[ProcNum - 1];

		for (int i = 1; i < ProcNum; i++)
			displs[i] = displs[i - 1] + sendCount[i - 1];

		
		for (int i = 0; i < ProcNum; i++)
		{
			recvCount[i] = sendCount[i] / n;
			gatherDispls[i] = displs[i] / n;
		}

	}

	MPI_Bcast(multVector,    n,       MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&border,       1,       MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(boolSendCount, 2,       MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs,        ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(recvCount,     ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(gatherDispls,  ProcNum, MPI_INT, 0, MPI_COMM_WORLD);

	int* bufferM = new int[boolSendCount[ProcRank >= border]];
	for (int i = 0; i < boolSendCount[ProcRank >= border]; i++)
		bufferM[i] = 0;

	MPI_Scatterv(matrix, sendCount, displs, MPI_INT, bufferM, boolSendCount[ProcRank >= border], MPI_INT, 0, MPI_COMM_WORLD);

	int* vectorElem = new int[recvCount[ProcRank]];
	for (int i = 0; i < recvCount[ProcRank]; i++)
		vectorElem[i] = 0;
	for (int i = 0; i < recvCount[ProcRank]; i++)
		for (int j = 0; j < n; j++)
			vectorElem[i] += bufferM[i * n + j] * multVector[j];

	MPI_Gatherv(vectorElem, recvCount[ProcRank], MPI_INT, resultVectorP, recvCount, gatherDispls, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		cout << "Result vector " << "                      " << "Difference vector" << endl;
		for (int i = 0; i < n; i++)
			cout << resultVectorP[i] << "                                " << resultVectorP[i] - resultVectorC[i] << endl;
		
		delete[] matrix;
		delete[] resultVectorP;
		delete[] resultVectorC;
		delete[] sendCount;
	}
	delete[] bufferM;
	delete[] vectorElem;
	delete[] multVector;
	delete[] displs;
	delete[] recvCount;
	delete[] gatherDispls;
	

	MPI_Finalize();
	return 0;
}