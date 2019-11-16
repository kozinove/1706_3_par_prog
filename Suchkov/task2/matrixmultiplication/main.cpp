#include <mpi.h>
#include <iostream>
#include <ctime>
#include <random>


using namespace std;

inline int* CreateMatrix(const int n, const int m)
{
	int* Matrix = new int[n * m];
	std::random_device rd;
	std::mt19937 mersenne(rd());

	for (int i = 0; i < n * m; i++)
	{
		Matrix[i] = (mersenne() % 10);
	}

	return Matrix;
}

inline void ShowMatrix(int* Matrix, const int m, const int n)
{
	int k = 1;
	for (int i = 0; i < n * m; i++, k++)
	{
		cout << Matrix[i] << "\t";
		if (k == m)
		{
			cout << endl;
			k = 0;
		}
	}
	cout << endl;
}


int main(int argc, char* argv[])
{
	int proc_id = 0, proc_num = 0;
	int* Matrix_A = nullptr;
	int* Matrix_B = nullptr;
	int* Matrix_C = nullptr;
	int row_A = 0, column_A = 0, row_B = 0, column_B = 0, column_C = 0, row_C = 0;
	int* Matrix_B_T = nullptr;
	int* receve_A = nullptr;
	int* receve_B = nullptr;
	int* tmp_C_part = nullptr;
	int how_to_count_c = 0;

	double time_start = 0;
	double time_work = 0;

	int* displs_A = nullptr;
	int* sendcounts_A = nullptr;
	int* displs_B = nullptr;
	int* sendcounts_B = nullptr;
	int* gather_displs = nullptr;
	int* gather_receve = nullptr;

	if (argc != 4)
	{
		cout << "\n\nERROR MUST BE 4 ARGUMETNS\n\n";
	}

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	/*try {
		column_A = atoi(argv[2]);
		row_C = row_A = atoi(argv[1]);
		column_C = column_B = atoi(argv[4]);
		row_B = atoi(argv[3]);
		if (row_B != column_A) row_B = column_A;
	}
	catch (...)
	{
		cout << "INVALID ARGUMENTS";
	}*/

	column_A = 3; // строка
	row_A = 3; // столбец
	column_B = 3; // строка
	row_B = 3; // столбец
	if (row_B != column_A) row_B = column_A;
	row_C = row_A;
	column_C = column_B;

	column_A = column_B = column_C = row_A = row_B = row_C = row_A; // debug


	displs_A = new int[proc_num];
	sendcounts_A = new int[proc_num];
	displs_B = new int[proc_num];
	sendcounts_B = new int[proc_num];
	gather_displs = new int[proc_num];
	gather_receve = new int[proc_num];

	for (int i = 0; i < proc_num; i++)
	{
		displs_A[i] = 0;
		sendcounts_A[i] = 0;
		displs_B[0] = 0;
		sendcounts_B[0] = 0;
	}

	if (proc_id == 0)
	{
		cout << "--------------------------------\n\n";
		Matrix_A = CreateMatrix(column_A, row_A);
		Matrix_B = CreateMatrix(column_B, row_B);
		cout << "MATRIX A \n";
		ShowMatrix(Matrix_A, column_A, row_A);
		cout << "MATRIX B \n";
		ShowMatrix(Matrix_B, column_B, row_B);
		Matrix_C = new int[column_C * row_C];
		for (int i = 0; i < column_C * row_C; i++)
			Matrix_C[i] = 0;
		cout << "--------------------------------\n\n";

		//-----------------------------------------------
		time_start = clock();

		for (int i = 0; i < row_A; i++)
			for (int j = 0; j < column_B; j++)
				for (int r = 0; r < row_B; r++)
					Matrix_C[i * column_C + j] += Matrix_A[i * column_A + r] * Matrix_B[r * row_B + j];

		time_work = clock() - time_start;
		cout << "\nMatrix C is = \n" << endl;
		ShowMatrix(Matrix_C, column_C, row_C);
		cout << "\nTime of sequential algorithm = " << time_work << " ms" << endl << endl;
		for (int i = 0; i < column_C * row_C; i++)
			Matrix_C[i] = 0;

		//-----------------------------------------------

		time_start = clock();
		Matrix_B_T = new int[row_C * column_C];
		for (int i = 0; i < row_B; i++)
			for (int j = 0; j < column_B; j++)
				Matrix_B_T[i * column_B + j] = Matrix_B[j * column_B + i];
		
		how_to_count_c = row_C * column_C / proc_num; // элементов матрицы C считает proc_id

		if (proc_num > row_C * column_C)
		{
			how_to_count_c = 1;
			sendcounts_A[0] = column_A * how_to_count_c / row_C ; //(how_to_count_c + row_C * column_C % proc_num) * column_A;
			displs_A[0] = 0;
			sendcounts_B[0] = row_B * how_to_count_c / column_C ; //(how_to_count_c + row_C * column_C % proc_num) * row_B;
			displs_B[0] = 0;
			gather_receve[0] = how_to_count_c;
			gather_displs[0] = 0;
		}
		else
		{
			sendcounts_A[0] = column_A * how_to_count_c / row_C + column_A * (row_C * column_C % proc_num) / row_C; //(how_to_count_c + row_C * column_C % proc_num) * column_A;
			displs_A[0] = 0;
			sendcounts_B[0] = row_B * how_to_count_c / column_C + row_B * (row_C * column_C % proc_num) / column_C; //(how_to_count_c + row_C * column_C % proc_num) * row_B;
			displs_B[0] = 0;
			gather_receve[0] = how_to_count_c + row_C * column_C % proc_num;
			gather_displs[0] = 0;
		}

		for (int i = 1; i < proc_num; i++)
		{
			gather_receve[i] = how_to_count_c;
			gather_displs[i] = gather_displs[i - 1] + gather_receve[i - 1];
		}

		for (int i = 1; i < proc_num; i++)
		{

			sendcounts_A[i] = column_A * how_to_count_c / row_C;   //how_to_count_c * column_A;
			displs_A[i] = sendcounts_A[i - 1] + displs_A[i - 1];
			sendcounts_B[i] = row_B * how_to_count_c / column_C;   //how_to_count_c * row_B;
			displs_B[i] = displs_B[i - 1] + displs_B[i - 1];

		}

	}
	
	MPI_Bcast(sendcounts_A, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs_A, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(sendcounts_B, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs_B, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(gather_receve, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(gather_displs, proc_num, MPI_INT, 0, MPI_COMM_WORLD);

	receve_A = new int[sendcounts_A[proc_id]];
	receve_B = new int[sendcounts_B[proc_id]];

	for (int i = 0; i < sendcounts_A[proc_id]; i++)	receve_A[i] = 0;
	for (int i = 0; i < sendcounts_B[proc_id]; i++) receve_B[i] = 0;

	MPI_Scatterv(Matrix_A, sendcounts_A, displs_A, MPI_INT, receve_A, sendcounts_A[proc_id], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(Matrix_B_T, sendcounts_B, displs_B, MPI_INT, receve_B, sendcounts_B[proc_id], MPI_INT, 0, MPI_COMM_WORLD);

	if (proc_id < row_C * column_C)
	{
		tmp_C_part = new int[how_to_count_c];
		for (int i = 0; i < how_to_count_c; i++)
			tmp_C_part[i] = i;



		/*for (int i = 0; i < how_to_count_c / column_C; i++)
			for (int j = 0; j < sendcounts_A[proc_id] / row_A; j++)
				for (int r = 0; r < sendcounts_B[proc_id] / column_B; r++)
					tmp_C_part[i * how_to_count_c + j] += receve_A[i * sendcounts_A[proc_id] / row_A + r] * receve_B[r * sendcounts_B[proc_id] / column_B + j];*/



		for (int i = 0; i < sendcounts_A[proc_id]; i++)
		{
			cout << receve_A[i] << "  ";
			cout << "proc_id = " << proc_id << "  ";
		}
		cout << endl << endl << endl;

		for (int i = 0; i < sendcounts_B[proc_id]; i++)
		{
			cout << receve_B[i] << "  ";
			cout << "proc_id = " << proc_id << "  ";
		}
		cout << endl << endl << endl << " ITOG\n";

		for (int i = 0; i < how_to_count_c; i++)
		{
			cout << tmp_C_part[i] << "  ";
			cout << "proc_id = " << proc_id << "  ";
		}
		cout << endl;
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Gatherv(tmp_C_part, how_to_count_c, MPI_INT, Matrix_C, gather_receve, gather_displs, MPI_INT, 0, MPI_COMM_WORLD);

	}

	if (proc_id == 0)
	{

		time_work = clock() - time_start;

		cout << endl << "\nMatrix is =\n";
		ShowMatrix(Matrix_C, column_C, row_C);
		cout << "Time of parallel algorithm = " << time_work << " ms" << endl << endl;
		cout << "--------------------------------\n\n";
	}

	delete[] receve_A;
	delete[] receve_B;
	delete[] sendcounts_A;
	delete[] displs_A;
	delete[] sendcounts_B;
	delete[] displs_B;
	delete[] Matrix_A;
	delete[] Matrix_B;
	delete[] Matrix_C;
	delete[] Matrix_B_T;
	delete[] tmp_C_part;
	delete[] gather_displs;
	delete[] gather_receve;

	return MPI_Finalize();;
}
