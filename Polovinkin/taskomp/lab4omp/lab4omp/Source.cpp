#include<omp.h>
#include <iostream>
#include <time.h>

using namespace std;
//создаем матрицы
double* CreateMatrix(int N)
{
	double* matrix = new double[N * N];

	return matrix;
}
//печатаем матрицу
void PrintMatrix(double* matrix, int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			cout << matrix[i * N + j] << " ";
		cout << endl;
	}

	cout << endl;
}
// заполняем матрицу
void GenerateRandomMatrix(double* matrix1, double* matrix2, int N)
{


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			matrix1[i * N + j] = (rand() % 101 - 50) / 10.0;
			matrix2[i * N + j] = (rand() % 101 - 50) / 10.0;
		}
	}
}
//делаем матрицу плотной
void Null(double* matrix1, double* matrix2, int N)
{
	for (int i = 0; i < 1; i++)
	{
		matrix1[i] = 0;
		matrix2[i] = 0;
	}
}
// обычный метод строчка на столбец
double* Simple_Mmult(double* matrix1, double* matrix2, int N) {
	double* Rez = CreateMatrix(N);
	int sum;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			sum = 0;
			for (int k = 0; k < N; k++)
				sum += matrix1[i * N + k] * matrix2[k * N + j];
			Rez[i * N + j] = sum;
		}
	return Rez;
}

// сумма 2-х матриц
double* Add(double* matrix1, double* matrix2, int N) {
	double* Rez = CreateMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j];

	return Rez;
}

// сумма 4-х матриц
double* Add(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N) {
	double* Rez = CreateMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j] + matrix3[i * N + j] + matrix4[i * N + j];

	return Rez;
}

// вычитание 2-ух матриц
double* Sub(double* matrix1, double* matrix2, int N) {
	double* Rez = CreateMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i * N + j] = matrix1[i * N + j] - matrix2[i * N + j];

	return Rez;
}

// сумма 3-х матриц вычесть 4-ую
double* Sub(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N)
{
	double* Rez = CreateMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j] + matrix3[i * N + j] - matrix4[i * N + j];

	return Rez;
}
// алгоритм Штрассена
double* Str_alg(double* matrix1, double* matrix2, int N, int threshold)
{
	double* Rez;

	if (N <= threshold) // меньше критического значения
		Rez = Simple_Mmult(matrix1, matrix2, N);
	else
	{
		Rez = CreateMatrix(N);
		N = N / 2;// разбили матрицу на 4 размеров Н/2 каждая

		double* A[4]; double* B[4]; double* C[4]; double* P[7];

		double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
		double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

		// выделение памяти под вспомогательные матрицы 
		for (int i = 0; i < 4; i++)
		{
			A[i] = CreateMatrix(N);
			B[i] = CreateMatrix(N);
		}

		// разбиение матрицы на 4 блока 
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				A[0][i * N + j] = matrix1[2 * i * N + j];
				A[1][i * N + j] = matrix1[2 * i * N + j + N];
				A[2][i * N + j] = matrix1[2 * i * N + j + 2 * N * N];
				A[3][i * N + j] = matrix1[2 * i * N + j + 2 * N * N + N];

				B[0][i * N + j] = matrix2[2 * i * N + j];
				B[1][i * N + j] = matrix2[2 * i * N + j + N];
				B[2][i * N + j] = matrix2[2 * i * N + j + 2 * N * N];
				B[3][i * N + j] = matrix2[2 * i * N + j + 2 * N * N + N];
			}

		// выполнение рекурсивного умножения 
		TMP1 = Add(A[0], A[3], N);
		TMP2 = Add(B[0], B[3], N);
		P[0] = Str_alg(TMP1, TMP2, N, threshold); // (A11 + A22)*(B11 + B22)

		TMP3 = Add(A[2], A[3], N);
		P[1] = Str_alg(TMP3, B[0], N, threshold); // (A21 + A22)*B11

		TMP4 = Sub(B[1], B[3], N);
		P[2] = Str_alg(A[0], TMP4, N, threshold); // A11*(B12 - B22)

		TMP5 = Sub(B[2], B[0], N);
		P[3] = Str_alg(A[3], TMP5, N, threshold); // A22*(B21 - B11)

		TMP6 = Add(A[0], A[1], N);
		P[4] = Str_alg(TMP6, B[3], N, threshold); // (A11 + A12)*B22

		TMP7 = Sub(A[2], A[0], N);
		TMP8 = Add(B[0], B[1], N);
		P[5] = Str_alg(TMP7, TMP8, N, threshold); // (A21 - A11)*(B11 + B12)

		TMP9 = Sub(A[1], A[3], N);
		TMP10 = Add(B[2], B[3], N);
		P[6] = Str_alg(TMP9, TMP10, N, threshold); // (A12 - A22)*(B21 + B22)

		// нахождение результирующих значений блоков 
		C[0] = Sub(P[0], P[3], P[6], P[4], N); // P1 + P4 - P5 + P7
		C[1] = Add(P[2], P[4], N); // P3 + P5
		C[2] = Add(P[1], P[3], N); // P2 + P4
		C[3] = Sub(P[0], P[2], P[5], P[1], N); // P1 - P2 + P3 + P6

		// формирование результирующей матрицы 
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				Rez[i * 2 * N + j] = C[0][i * N + j];
				Rez[i * 2 * N + j + N] = C[1][i * N + j];
				Rez[i * 2 * N + j + 2 * N * N] = C[2][i * N + j];
				Rez[i * 2 * N + j + 2 * N * N + N] = C[3][i * N + j];
			}

		// освобождение памяти 
		for (int i = 0; i < 4; i++) {
			delete[] A[i];
			delete[] B[i];
			delete[] C[i];
		}

		for (int i = 0; i < 7; i++) {
			delete[] P[i];
		}

		delete[]TMP1; delete[]TMP2; delete[]TMP3; delete[]TMP4; delete[]TMP5;
		delete[]TMP6; delete[]TMP7; delete[]TMP8; delete[]TMP9; delete[]TMP10;
	}

	return Rez;
}
// OMP алгоритм Штрассена
double* Shtr_alg_omp(double* matrix1, double* matrix2, int N, int threshold)
{
	double* Rez;

	if (N <= threshold) {
		Rez = Simple_Mmult(matrix1, matrix2, N);
	}
	else {
		Rez = CreateMatrix(N);
		N = N / 2;
		double* A[4];
		double* B[4];
		double* C[4];
		double* P[7];

		for (int i = 0; i < 4; i++) {
			A[i] = CreateMatrix(N);
			B[i] = CreateMatrix(N);
		}
#pragma omp parallel
		{
			int i, j;
// распараллеливаем цикл for
#pragma omp for private(i,j) 
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)
				{
					int index_new = i * N + j, index_old = 2 * i * N + j, N_N = 2 * N * N;
					A[0][index_new] = matrix1[index_old];
					A[1][index_new] = matrix1[index_old + N];
					A[2][index_new] = matrix1[index_old + N_N];
					A[3][index_new] = matrix1[index_old + N_N + N];

					B[0][index_new] = matrix2[index_old];
					B[1][index_new] = matrix2[index_old + N];
					B[2][index_new] = matrix2[index_old + N_N];
					B[3][index_new] = matrix2[index_old + N_N + N];
				}
			// одна секция выполняется одним потоком
#pragma omp sections 
			{
double* TMP = Add(A[0], A[3], N);
					double* _TMP = Add(B[0], B[3], N);
#pragma omp section
				{
					
					P[0] = Simple_Mmult(TMP, _TMP, N);
					delete[] TMP;
					delete[] _TMP;
				}
#pragma omp section
				{
					TMP = Add(A[2], A[3], N);
					P[1] = Simple_Mmult(TMP, B[0], N);
					delete[] TMP;
				}
#pragma omp section
				{
					TMP = Sub(B[1], B[3], N);
					P[2] = Simple_Mmult(A[0], TMP, N);
					delete[] TMP;
				}
#pragma omp section
				{
					TMP = Sub(B[2], B[0], N);
					P[3] = Simple_Mmult(A[3], TMP, N);
					delete[] TMP;
				}
#pragma omp section
				{
					TMP = Add(A[0], A[1], N);
					P[4] = Simple_Mmult(TMP, B[3], N);
					delete[] TMP;
				}
#pragma omp section
				{
					TMP = Sub(A[2], A[0], N);
					_TMP = Add(B[0], B[1], N);
					P[5] = Simple_Mmult(TMP, _TMP, N);
					delete[] TMP;
					delete[] _TMP;
				}
#pragma omp section
				{
					TMP = Sub(A[1], A[3], N);
					_TMP = Add(B[2], B[3], N);
					P[6] = Simple_Mmult(TMP, _TMP, N);
					delete[] TMP;
					delete[] _TMP;
				}
			}
			// одна секция выполняется одним потоком
#pragma omp sections
			{
#pragma omp section
				C[0] = Sub(P[0], P[3], P[6], P[4], N);

#pragma omp section
				C[1] = Add(P[2], P[4], N);

#pragma omp section
				C[2] = Add(P[1], P[3], N);

#pragma omp section
				C[3] = Sub(P[0], P[2], P[5], P[1], N);
			}
// распараллеливаем цикл
#pragma omp for private(i,j)
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++) {
					Rez[i * 2 * N + j] = C[0][i * N + j];
					Rez[i * 2 * N + j + N] = C[1][i * N + j];
					Rez[i * 2 * N + j + 2 * N * N] = C[2][i * N + j];
					Rez[i * 2 * N + j + 2 * N * N + N] = C[3][i * N + j];
				}
		}

		for (int i = 0; i < 4; i++) {
			delete[] A[i];
			delete[] B[i];
			delete[] C[i];
		}

		for (int i = 0; i < 7; i++) {
			delete[] P[i];
		}


	}

	return Rez;
}
int main()
{
	double* matr_A = NULL;
	double* matr_B = NULL;
	double* matr_Rez_Str = NULL;
	double* matr_Rez_Str2 = NULL;
	double* matr_Rez_Omp = NULL;
	double StartStrAlg = 0;	// время старта последовательного алгоритма Штрассена
	double TimeStrAlg = 0;// время, затраченное на работу последовательной версии алгоритма Штрассена
	int thr = 64;
	int k = 9;
	double StartStrAlg2 = 0;// время старта последовательного обычного алгоритма 
	double TimeStrAlg2 = 0;// время, затраченное на работу обычной версии алгоритма
	double EndOmpAlg = 0;// Время окончания omp алгоритма
	double StartOmpAlg = 0;// Время старта omp алгоритма
	double TimeOmpAlg = 0;// Время,затраченно на работу omp алгоритма



	int N = (int)pow(2.0, k);
	matr_A = CreateMatrix(N);
	matr_B = CreateMatrix(N);

	GenerateRandomMatrix(matr_A, matr_B, N);
	Null(matr_A, matr_B, N);

// обычный алгоритм
	StartStrAlg2 = clock();

	matr_Rez_Str2 = Simple_Mmult(matr_A, matr_B, N);

	TimeStrAlg2 = clock() - StartStrAlg2;

	cout << "Usually linear:" << TimeStrAlg2 << endl;

	// последовательный алгоритм 
	StartStrAlg = clock();

	matr_Rez_Str = Str_alg(matr_A, matr_B, N, thr);

	TimeStrAlg = clock() - StartStrAlg;

	cout << "Strassen linear:" << TimeStrAlg << endl;

	// omp алгоритм
	StartOmpAlg = omp_get_wtime();

	matr_Rez_Omp = Shtr_alg_omp(matr_A, matr_B, N, thr);

	EndOmpAlg = omp_get_wtime();

	TimeOmpAlg = EndOmpAlg - StartOmpAlg;

	cout << "Strassen OpenMp:" << TimeOmpAlg << endl;

	delete[] matr_Rez_Omp;
	delete[] matr_Rez_Str;
	delete[] matr_Rez_Str2;
	delete[] matr_B;
	delete[] matr_A;

	return 0;

}