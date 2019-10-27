#include <mpi.h>
#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
	int ProcRank = 0, ProcNum = 0;
	int* Matrix = nullptr;
	int* Matrix_vector = nullptr;
	int* Get_vector = nullptr;
	int Receive_vector_size = 0;
	int Partial_max = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int Columns = ProcNum, Rows = ProcNum;
	int* Max = new int[Rows];
	srand(unsigned(time(0)));

	if (ProcRank == 0) {
		Matrix = new int[Columns * Rows];
		Matrix_vector = new int[Columns * Rows];
		for (int i = 0; i < Columns * Rows; ++i) { //matrix init 

			Matrix[i] = (rand() % 1000);

		}
		for (int i = 0; i < Rows; i++)
		{
			for (int j = 0; j < Columns; j++)
			{
				cout << Matrix[i * Columns + j] << " ";
				cout << " ";
			}
			cout << endl;
		}
		cout << endl;

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
	Receive_vector_size = Columns;
	int* Displs = new int[ProcNum]; //mpi 
	int* Sendcounts = new int[ProcNum]; //mpi

	for (int i = 0; i < ProcNum; i++) {
		Sendcounts[i] = Receive_vector_size;
		Displs[i] = (i * Receive_vector_size);
	}

	MPI_Bcast(&Receive_vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

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

	MPI_Gather(&Partial_max_part, 1, MPI_INT, &Get_vector[ProcRank], 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		for (int i = 0; i < Rows; i++) {
			cout << "\nMax of " << i + 1 << " row in matrix = " << Max[i] << endl;
			if (Max[i] > Partial_max)
				Partial_max = Max[i];
		}
		cout << endl;
		cout << "Max of all elements = " << Partial_max << endl;
		delete[] Max;
		delete[] Matrix;
		delete[] Matrix_vector;
		delete[] Sendcounts;
		delete[] Displs;
	}

	delete[] Get_vector;
	MPI_Finalize();
	return 0;
}