#include <mpi.h>
#include <iostream>
#include <ctime>
#include <vector>
#include <windows.h>


using namespace std;

inline int** CreateMatrix(const int n, const int m)
{
	int** Matrix = new int*[n];
	for (int i = 0; i < n; ++i)
		Matrix[i] = new int[m];

	srand(unsigned(time(0)));

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			Matrix[i][j] = (rand() % 100) - 50; // îò -49 äî 49
	}

	return Matrix;
}

inline void ShowMatrix(int** Matrix, const int n, const int m)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			cout << Matrix[i][j] << " ";

		cout << endl;
	}
}

int main(int argc, char* argv[])
{
	int proc_id = 0, proc_num = 0;
	int** Matrix = nullptr;
	int* Matrix_vector = nullptr;
	int m, n;
	int* receve_vector = nullptr;
	int work_size = 0;

	double time_start = 0;
	double parall_time_work = 0;
	double sequent_time_work = 0;
	int parallel_sum = 0;
	int sequential_sum = 0;
	int* displs = nullptr;
	int* sendcounts = nullptr;
	int partial_sum_part = 0;
	MPI_Status Status;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	m = atoi(argv[2]);
	n = atoi(argv[1]);
	
	displs = new int[proc_num];
	sendcounts = new int[proc_num];

	if (proc_id == 0)
	{
		Matrix = CreateMatrix(n, m);
		ShowMatrix(Matrix, n, m);
		
		//-----------------------------------------------
		time_start = clock();
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				sequential_sum += Matrix[i][j];
		}
		sequent_time_work = clock() - time_start;
		cout << "\nSum of all elements in matrix is = " << sequential_sum << endl;
		cout << "\nTime of sequential algorithm = " << sequent_time_work << " ms" << endl << endl;

		//-----------------------------------------------

		time_start = clock();

		Matrix_vector = new int[n * m];
		for (int i = 0, k = 0; i < n; i++)
			for (int j = 0; j < m; j++, k++)
				Matrix_vector[k] = Matrix[i][j];

		sendcounts[0] = (n * m) / proc_num + (n * m) % proc_num;
		displs[0] = 0;
		for (int i = 1; i < proc_num; i++)
		{
			sendcounts[i] = (n * m) / proc_num;
			displs[i] = sendcounts[i - 1] + displs[i - 1];
		}
	}
	
	MPI_Bcast(sendcounts, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, proc_num, MPI_INT, 0, MPI_COMM_WORLD);

	work_size = sendcounts[proc_id];

	receve_vector = new int[work_size];
	for (int i = 0; i < work_size; i++)
	{
		receve_vector[i] = 0;
	}

	MPI_Scatterv(Matrix_vector, sendcounts, displs, MPI_INT, receve_vector, work_size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < work_size; i++)
	{
		cout << receve_vector[i] << "  ";
	}
	cout << endl;

	if (work_size * proc_id < m * n)
	{
		for (int i = 0; i < work_size; i++)
		{
			partial_sum_part += receve_vector[i];
			cout << proc_id << "   " << partial_sum_part << "  " << receve_vector[i] << endl;
		}
		cout << "partial_sum_part[" << proc_id << "] = " << partial_sum_part << endl;
	}
	MPI_Reduce(&partial_sum_part, &parallel_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (proc_id == 0)
	{
		
		parall_time_work = clock() - time_start;

		cout << endl << "Sum of all elements in matrix is = " << parallel_sum << endl;
		cout << "Time of parallel algorithm = " << parall_time_work << " ms" << endl << endl;


		for (int i = 0; i < n; i++)
			delete[] Matrix[i];
		delete[] Matrix_vector;
		//delete[] sendcounts;
		//delete[] displs;
		

	}
	delete[] receve_vector;
	//delete[] partial_sum_part;
	delete[] sendcounts;
	delete[] displs;

	return MPI_Finalize();;
}
