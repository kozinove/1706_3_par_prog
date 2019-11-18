#include <stdio.h>
#include <mpi.h>
#include <math.h>
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

	if (argc < 5)
	{
		if(ProcRank == 0)
			cout << "NO ENOUGH ARGUMENTS, WILL BE USED SERVICE DATA" << endl;
	}
	if (argc > 5)
	{
		if (ProcRank == 0)
			cout << "TOO MUCH ARGUMENTS, SOME WON`T BE USED" << endl;
	}
	
	int* matrix        = nullptr;
	int* sendCount     = nullptr;
	int* resultVectorP = nullptr;  //Parallel   mult
	int* resultVectorC = nullptr;  //Consistent mult
	int Cols = 5;
	if (argc > 1)
	{ 
		Cols = atoi(argv[1]);
	}
	int Rows = Cols;
	if (argc > 2)
	{
		Rows = atoi(argv[2]);
	}

	if (Cols < 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED COLS = |COLS|" << endl;
		Cols = abs(Cols);
	}
	else if (Cols == 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED COLS = 5" << endl;
		Cols = 5;
	}
	if (Rows < 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED ROWS = |ROWS|" << endl;
		Rows = abs(Rows);
	}
	else if (Rows == 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED ROWS = COLS" << endl;
		Rows = Cols;
	}
	int* multVector   = new int[Cols];
	int* displs       = new int[ProcNum];
	int* recvCount    = new int[ProcNum];
	int* gatherDispls = new int[ProcNum];
	int border = 0;
	int boolSendCount[2];
	if (ProcRank == 0)
	{
		int Min = 0;
		int Max = 9;
		if (argc > 3)
			Min  = atoi(argv[3]);
		if (argc > 5)
			Max  = atoi(argv[4]);
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
		cout << "USING DATA: " << endl << "Cols = " << Cols << endl << "Rows = " << Rows << endl << "Min = " << Min << endl << "Max = " << Max << endl;

	
		srand(time(NULL));
		matrix = new int[Cols * Rows];

		for (int i = 0; i < Cols; i++)
		{
			multVector[i] = Min + rand() % (Max - Min + 1);
			for (int j = 0; j < Rows; j++)
				matrix[i + j*Cols] = Min + rand() % (Max - Min + 1);
		}	

		cout << "Matrix: " << endl;
		for (int i = 0; i < Rows; i++)
		{
			for (int j = 0; j < Cols; j++)
				cout << matrix[i * Cols + j] << "   ";
			cout << endl;
		}
		cout << "Vector: " << endl;
		for (int i = 0; i < Cols; i++)
			cout << matrix[i] << endl;

		cout << endl << "--------------------------------" << endl << endl;
			
		resultVectorC = new int[Rows];
		for (int i = 0; i < Rows; i++)
			resultVectorC[i] = 0;
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				resultVectorC[i] += matrix[i * Cols + j] * multVector[j];      
			
		resultVectorP = new int[Rows];
		for (int i = 0; i < Rows; i++)
			resultVectorP[i] = 0;
				
		int startNum = Rows / ProcNum;
		sendCount = new int[ProcNum];
		for (int i = 0; i < ProcNum; i++)
		{
			sendCount[i] = startNum * Cols; ////
			displs[i] = 0;
		}
		border = Rows - startNum * ProcNum;
		for (int i = 0; i < border; i++)
			sendCount[i] += Cols; ////

		boolSendCount[0] = sendCount[0];
		boolSendCount[1] = sendCount[ProcNum - 1];

		for (int i = 1; i < ProcNum; i++)
			displs[i] = displs[i - 1] + sendCount[i - 1];

		
		for (int i = 0; i < ProcNum; i++)
		{
			recvCount[i] = sendCount[i] / Cols;
			gatherDispls[i] = displs[i] / Cols;
		}

	}

	MPI_Bcast(multVector,    Cols,    MPI_INT, 0, MPI_COMM_WORLD);
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
		for (int j = 0; j < Cols; j++)
			vectorElem[i] += bufferM[i * Cols + j] * multVector[j];


	MPI_Gatherv(vectorElem, recvCount[ProcRank], MPI_INT, resultVectorP, recvCount, gatherDispls, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		cout << "Result vector " << "                      " << "Difference vector" << endl;
		for (int i = 0; i < Rows; i++)
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
