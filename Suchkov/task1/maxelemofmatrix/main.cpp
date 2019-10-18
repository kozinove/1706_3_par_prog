#include <mpi.h>
#include <iostream>
#include <ctime>
#include <vector>
#include <windows.h>


using namespace std;

inline int** CreateMatrix(const int n, const int m)
{
	int** Matrix = new int* [n];
	for (int i = 0; i < n; ++i)
		Matrix[i] = new int[m];

	srand(unsigned(time(0)));

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			Matrix[i][j] = (rand() % 100) - 50; // от -49 до 49
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
	int* partial_sum = 0;
	MPI_Status Status;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	
	m = atoi(argv[2]);
	n = atoi(argv[1]);
	
	//cout << proc_num << "   " << proc_id << endl;
	if (proc_id == 0)
	{
		Matrix = CreateMatrix(n, m);
		ShowMatrix(Matrix, n, m);
		partial_sum = new int[proc_num];
		for (int i = 0; i < proc_num; i++)
			partial_sum[i] = 0;
		
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

	}
	work_size = (n * m) / proc_num;

	receve_vector = new int[work_size + 3];
	for (int i = 0; i < work_size + 3; i++)
		receve_vector[i] = 0;

	int* displs = new int[proc_num];
	/*for (int i = 0; i < proc_num; i++)
		displs[i] = (i * work_size);*/

	int* sendcounts = new int[proc_num];
	for (int i = 0; i < proc_num; i++)
	{
		sendcounts[i] = work_size;
		displs[i] = (i * work_size);
		if ((work_size % proc_num != 0) && (i == (proc_num - 1)))
		{
			sendcounts[i] = work_size + 1;
			displs[i]++;
		}
	}

	int* partial_sum_part = new int[work_size];
	for (int i = 0; i < work_size; i++)
		partial_sum_part[i] = 0;

	MPI_Scatterv(Matrix_vector, sendcounts, displs, MPI_INT, receve_vector, work_size + 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < work_size; i++)
	{
		partial_sum_part[proc_id] += receve_vector[i];
		cout << proc_id << "   " << partial_sum_part[proc_id] << "  " << receve_vector[i] << endl;
	}
	cout << "partial_sum_part[" << proc_id << "] = " << partial_sum_part[proc_id] << endl;

	MPI_Reduce(partial_sum_part, partial_sum, proc_num, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (proc_id == 0)
	{
		for (int i = 0; i < proc_num; i++)
			parallel_sum += partial_sum[i];
		parall_time_work = clock() - time_start;

		cout << endl << "Sum of all elements in matrix is = " << parallel_sum << endl;
		cout << "Time of parallel algorithm = " << parall_time_work << " ms" << endl << endl;
		
		
		for (int i = 0; i < n; i++)
				delete[] Matrix[i];
		delete[] Matrix_vector;
		delete[] receve_vector;
		delete[] partial_sum;
		delete[] partial_sum_part;
		delete[] sendcounts;
		delete[] displs;

	}

	//MPI_Finalize();
	
	return MPI_Finalize();;
}

