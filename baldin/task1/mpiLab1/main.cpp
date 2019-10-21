#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	int ProcNum = 0;
	int ProcRank = 0;
	int n = 0;
	int startNum = 0;
	int BufferMem = 0;
	int CMax = 0; //Consistent maximum 
	int PMax = 0; //Parallel maximum
	int Min = atoi(argv[2]);
	int Max = atoi(argv[3]);
	int* Matrix = nullptr;
	int* MaxArray = nullptr;
	int* displs = nullptr;
	int* sendCount = nullptr;
	

	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);	
	if (ProcRank == 0)
	{
		n = atoi(argv[1]);
		if (ProcNum > n * n)
		{
			Matrix = new int[n * n + ProcNum - n * n];
			for (int i = n * n; i < n * n + ProcNum - n * n; i++)
				Matrix[i] = -2147483647;
		}
		else
			Matrix = new int[n * n];
		srand(time(NULL));
		if (Max < Min)
		{
			int temp = Max;
			Max = Min;
			Min = temp;
		}
		for (int i = 0; i < n*n; i++)
				Matrix[i] = Min+rand()%(Max-Min);
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				cout << Matrix[i * n + j] << "   ";
			cout << endl;
		}
		cout << endl;

		MaxArray = new int[ProcNum];
		for (int i = 0; i < ProcNum; i++)
			MaxArray[i] = 0;

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (CMax < Matrix[i * n + j])
					CMax = Matrix[i * n + j];
		cout << "Consistent maximum = " << CMax << endl;

		sendCount = new int[ProcNum];
		startNum = n * n / ProcNum;
		for (int i = 0; i < ProcNum; i++)
			sendCount[i] = startNum;
		for (int i = 0; i < n * n - startNum * ProcNum; i++)
			sendCount[i]++;


		for (int i = 0; i < ProcNum; i++)
			if (sendCount[i] == 0)
				sendCount[i] = 1;

		displs = new int[ProcNum];
		for (int i = 0; i < ProcNum; i++)
			displs[i] = 0;
		for (int i = 1; i < ProcNum; i++)
			displs[i] = displs[i-1] + sendCount[i-1];

		BufferMem = sendCount[0];
	}
	MPI_Bcast(&BufferMem, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int* Buffer = new int[BufferMem];
	for (int i = 0; i < BufferMem; i++)
		Buffer[i] = -2147483647;

	MPI_Scatterv(Matrix, sendCount, displs, MPI_INT, Buffer, BufferMem, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 0; i < BufferMem; i++)
		if (PMax < Buffer[i])
			PMax = Buffer[i];
	MPI_Gather(&PMax, 1, MPI_INT, &MaxArray[ProcRank], 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		for (int i = 0; i < ProcNum; i++)
			if (PMax < MaxArray[i])
				PMax = MaxArray[i];
		cout << "Parallel max = " << PMax << "    PMax - CMax = " << PMax - CMax << endl;
		delete[] Matrix;
		delete[] MaxArray;
		delete[] sendCount;
		delete[] displs;

	}
	delete[] Buffer;

	MPI_Finalize();


	return 0;
}