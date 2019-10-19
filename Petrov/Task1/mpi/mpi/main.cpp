#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
	int ProcRank = 0, ProcNum = 0;
	int Columns = 10, Rows = 10;
	int* Matrix = nullptr;
	int* Matrix_vector = nullptr;
	int* Get_vector = nullptr;
	int Receive_vector_size = 0;
	int* Max = new int[Rows];
	int Partial_max = 0;
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		Matrix = new int[Columns * Rows];
		Matrix_vector = new int[Columns * Rows];
		for (int i = 0; i < Rows; ++i) { //matrix init 
			for (int j = 0; j < Columns; ++j)
			{
				Matrix[i * Columns + j] = (rand() % 100);
			}
		}
		for (int i = 0; i < Rows; i++) {
			Max[i] = 0;
			for (int j = 0; j < Rows; j++) {
				if (Matrix[i * Columns + j] > Max[i])
					Max[i] = Matrix[i * Columns + j];
			}
		}
		for (int i = 0, k = 0; i < Rows; i++)
			for (int j = 0; j < Columns; j++, k++)
				Matrix_vector[k] = Matrix[i * Columns + j];
	}
	Receive_vector_size  = (Rows * Columns) / ProcNum; 
	int* Displs = new int[ProcNum]; //mpi 
	int* Sendcounts = new int[ProcNum]; //mpi

	for (int i = 0; i < ProcNum; i++) {
		Sendcounts[i] = Receive_vector_size;
		Displs[i] = (i * Receive_vector_size);
	}
	Get_vector = new int[Receive_vector_size];
	for (int i = 0; i < Receive_vector_size; i++) {
		Get_vector[i] = 0;
	}
	MPI_Scatterv(Matrix_vector, Sendcounts, Displs, MPI_INT, Get_vector, Receive_vector_size, MPI_INT, 0, MPI_COMM_WORLD);

	int Partial_max_part = 0;

	for (int i = 0; i < Receive_vector_size; i++) {
		if (Get_vector[i] > Partial_max_part)
			Partial_max_part = Get_vector[i];
		cout << "ProcRank = " << ProcRank << " | Max elem = " << Partial_max_part << " | Current elem = " << Get_vector[i] << endl;
	}
	cout << "Partial_max_part[" << ProcRank << "] = " << Partial_max_part << endl;

	MPI_Reduce(&Partial_max_part, &Partial_max, ProcNum, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		for (int i = 0; i < Rows; i++)
			cout << "\nMax of "<< i+1 <<" row in matrix = " << Max[i] << endl;
		cout << endl;
		cout << "Partial_max: " << Partial_max << endl;
		delete[] Max;
		delete[] Matrix;
		delete[] Matrix_vector;
		delete[] Get_vector;
		delete[] Sendcounts;
		delete[] Displs;
	}
	MPI_Finalize();
	return 0;
}