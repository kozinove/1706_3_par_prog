#include <mpi.h>
#include <iostream>
#include <ctime>

using namespace std;

int process_count = 0; //Кол-во процессов
int cur_rank = 0;	   //Текущий процесс
int* A = 0, * B = 0, * C = 0, * A1 = 0, * B1 = 0, * C1 = 0, * res; //Все матрицы
int matrix_size = 0, block_size = 0;  //размер матрицы и блока
//MPI_Comm GridComm, ColComm, RowComm; //Новые коммуникаторы
MPI_Datatype MPI_BLOCK;  //Тип для пересылки блоков

int iter = 0; //итерация

//Инициализация матриц
void initialize()
{
	block_size = matrix_size / process_count; //вычисление размера блока
	//Выделение памяти на матрицы (полосы)
	A1 = new int[block_size * matrix_size];
	B1 = new int[block_size * matrix_size];
	if (cur_rank == 0)
	{
		// задание элементов исходных матриц
		srand(time(0));
		A = new int[matrix_size * matrix_size];
		B = new int[matrix_size * matrix_size];
		C = new int[matrix_size * matrix_size];
		res = new int[matrix_size * matrix_size];
		//случайная генерация исходных матриц
		for (int i = 0; i < matrix_size; ++i)
			for (int j = 0; j < matrix_size; ++j)
			{
				A[i * matrix_size + j] = rand() % 50;
				B[i * matrix_size + j] = rand() % 50;
				res[i * matrix_size + j] = 0;
			}
		C1 = C;
	}
	else
	{
		C1 = new int[block_size * matrix_size];
	}
	//Обнуление C1
	for (int i = 0; i < block_size * matrix_size; ++i)
		C1[i] = 0;
}
//Последовательный алгоритм
void serialMult()
{
	for (int i = 0; i < matrix_size; ++i)
		for (int k = 0; k < matrix_size; ++k)
			for (int j = 0; j < matrix_size; ++j)
				res[i * matrix_size + j] += A[i * matrix_size + k] * B[k * matrix_size + j];
}
//Базовая рассылка матриц
void matrixSending()
{
	MPI_Type_vector(matrix_size, block_size, matrix_size, MPI_INT, &MPI_BLOCK); //Блок матрицы
	MPI_Type_commit(&MPI_BLOCK);
	if (cur_rank == 0)
	{
		for (int r = 1; r < process_count; ++r)
		{
			//Отправка соответствующих блоков
			MPI_Send(A + r * matrix_size * block_size, matrix_size * block_size, MPI_INT, r, 0, MPI_COMM_WORLD);
			MPI_Send(B + r * block_size, 1, MPI_BLOCK, r, 0, MPI_COMM_WORLD);
		}
		for (int i = 0; i < block_size * matrix_size; ++i)
			A1[i] = A[i];
		for (int i = 0; i < matrix_size; ++i)
			for (int j = 0; j < block_size; ++j)
				B1[i * block_size + j] = B[i * matrix_size + j];
	}
	else
	{
		//Прием блоков
		MPI_Status s;
		MPI_Recv(A1, block_size * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
		MPI_Recv(B1, block_size * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
	}
}
//Перемножение блока B на  блок A
void multA1_B1()
{
	//Циклы в оптимальной последовательности
	for (int i = 0; i < block_size; ++i)
		for (int k = 0; k < matrix_size; k++)
			for (int j = 0; j < block_size; ++j)
				C1[i * matrix_size + j + block_size * ((iter + cur_rank) % process_count)] += A1[i * matrix_size + k] * B1[k * block_size + j];
}
//Рассылка B1
void B1_sending()
{
	MPI_Status s;
	int next_proc = (cur_rank + 1) % process_count;
	int prev_proc = (cur_rank + process_count - 1) % process_count;
	MPI_Sendrecv_replace(B1, block_size * matrix_size, MPI_INT, prev_proc, 0, next_proc, 0, MPI_COMM_WORLD, &s);
}
//Умножение
void mult()
{
	for (iter = 0; iter < process_count; ++iter)
	{
		multA1_B1();
		if (iter < process_count - 1)
			B1_sending();
	}
}
//Сбор результатов
void uniteResult()
{
	if (cur_rank == 0)
	{
		MPI_Status s;
		for (int r = 1; r < process_count; ++r)
		{
			MPI_Recv(C + r * matrix_size * block_size, matrix_size * block_size, MPI_INT, r, 0, MPI_COMM_WORLD, &s);
		}
	}
	else
		MPI_Send(C1, block_size * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
}
//Освобождение памяти
void finalize()
{
	delete[] A1;
	delete[] B1;
	delete[] C1;
	if (cur_rank == 0)
	{
		delete[] A;
		delete[] B;
		delete[] res;
	}
	MPI_Type_free(&MPI_BLOCK);
}

int main(int argc, char** argv)
{
	double start_time, end_time;
	if (argc > 1)
		matrix_size = atoi(argv[1]);
	else
		matrix_size = 5;
	//MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &cur_rank);

	initialize();
	if (cur_rank == 0)
	{
		start_time = MPI_Wtime();
		serialMult();
		end_time = MPI_Wtime();
		cout << "Serial: " << end_time - start_time << endl;
	}
	start_time = MPI_Wtime();
	matrixSending();
	mult();
	uniteResult();
	if (cur_rank == 0)
	{
		end_time = MPI_Wtime();
		cout << "Parallel: " << end_time - start_time << endl;
		bool flag = true;
		for (int i = 0; i < matrix_size * matrix_size; ++i)
			flag = flag && (C[i] == res[i]);
		if (flag)
			cout << "equal" << endl;
		else
			cout << "is not equal" << endl;
		//for (int i = 0; i < matrix_size; ++i)
		//{
		//	for (int j = 0; j < matrix_size; ++j)
		//		cout << C[i*matrix_size + j] << " ";
		//	cout << endl;
		//}
		//for (int i = 0; i < matrix_size; ++i)
		//{
		//	for (int j = 0; j < matrix_size; ++j)
		//		cout << res[i*matrix_size + j] << " ";
		//	cout << endl;
		//}
	}
	finalize();

	MPI_Finalize();
	return 0;
}