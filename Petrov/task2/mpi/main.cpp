#include <mpi.h>
#include <iostream>
#include <ctime>
using namespace std;

void PrintMatrix(int* Matrix, int Size)
{
	for (int i = 0; i < Size; ++i) {
		for (int j = 0; j < Size; ++j) {
			cout << Matrix[i * Size + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void InitMatrix(int* Matrix, int Rows, int Columns)
{
	for (int i = 0; i < Rows; ++i)
		for (int j = 0; j < Columns; ++j)
			Matrix[i * Columns + j] = (rand() % 10);
}

void Transpose(int* Matrix, int size)
{
	int t;
	for (int i = 0; i < size; ++i)
	{
		for (int j = i; j < size; ++j)
		{
			t = Matrix[i * size + j];
			Matrix[i * size + j] = Matrix[j * size + i];
			Matrix[j * size + i] = t;
		}
	}
}

int main(int argc, char* argv[]) {
	int Index;
	int NumOfProcess;
	int i, j, k, p;
	int ProcNum = 0;
	int ProcRank = 0;
	int Rows = 3;
	int Columns = 3;
	int Size = Rows;
	int* MatrixA = new int[Rows * Columns];
	int* MatrixB = new int[Rows * Columns];
	int* MatrixC = new int[Rows * Columns];
	int Temp = 0;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	srand(unsigned(time(0)));
	int ProcPart_A = Size / ProcNum;
	int ProcPart_B = Size / ProcNum;
	int Part_A = ProcPart_A * Size;
	int Part_B = ProcPart_B * Size;
	int* BufferA = new int[Part_A];
	int* BufferB = new int[Part_B];
	int* BufferC = new int[ProcPart_A * Size];
	if (ProcRank == 0) {
	for (i = 0; i < ProcPart_A * Size; i++)
		BufferC[i] = 0;
		InitMatrix(MatrixA, Rows, Columns);
		InitMatrix(MatrixB, Rows, Columns);
		PrintMatrix(MatrixA, Size);
		PrintMatrix(MatrixB, Size);
		Transpose(MatrixB, Size);
	}

	MPI_Scatter(MatrixA, Part_A, MPI_INT, BufferA, Part_A, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(MatrixB, Part_B, MPI_INT, BufferB, Part_B, MPI_INT, 0, MPI_COMM_WORLD);

	Temp = 0;
	for (i = 0; i < ProcPart_A; i++)
		for (j = 0; j < ProcPart_B; j++)
		{
			Temp = 0;
			for (k = 0; k < Size; k++)
				Temp += BufferA[i * Size + k] * BufferB[j * Size + k];
			BufferC[i * Size + j + ProcPart_A * ProcRank] = Temp;
		}
	int NextProc;
	int PrevProc;
	NextProc = ProcRank + 1;
	if (ProcRank == ProcNum - 1) NextProc = 0;
	PrevProc = ProcRank - 1;
	if (ProcRank == 0) PrevProc = ProcNum - 1;
	for (NumOfProcess = 0; NumOfProcess < ProcNum; NumOfProcess++)
	{
		MPI_Sendrecv_replace(BufferB, Part_B, MPI_INT, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);   // передача-принятие единого типа сообщений (рассылка B)
		Temp = 0;
		for (i = 0; i < ProcPart_A; i++)
		{
			for (j = 0; j < ProcPart_B; j++)
			{
				Temp = 0;
				for (k = 0; k < Size; k++)
					Temp += BufferA[i * Size + k] * BufferB[j * Size + k];
				if (ProcRank - NumOfProcess >= 0)
					Index = ProcRank - NumOfProcess;
				else
					Index = (ProcNum - NumOfProcess + ProcRank);
				BufferC[i * Size + j + Index * ProcPart_A] = Temp;
				 cout<< "ProcRank = " <<ProcRank<< " | Current elem = " << BufferC[i * Size + j + Index * ProcPart_A] << endl;
			}
		}
	}
	
	MPI_Gather(BufferC, ProcPart_A * Size, MPI_INT, MatrixC, ProcPart_A * Size, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		PrintMatrix(MatrixC, Size);
		delete[]MatrixA;
		delete[]MatrixB;
		delete[]MatrixC;
		delete[]BufferA;
		delete[]BufferB;
		delete[]BufferC;
	}
	MPI_Finalize();

	return 0;
}

//mpiexec -n 9 "C:\Users\\\\\\mpi.exe"