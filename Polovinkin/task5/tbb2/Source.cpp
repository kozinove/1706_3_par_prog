#include <iostream>
#include<math.h>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"
#include "tbb/tick_count.h"


using namespace std;
using namespace tbb;
//создаем матрицы
double* CreateMatrix(int N)
{
	double* matrix = new double[N * N];

	return matrix;
}

// заполн€ем матрицу
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
// алгоритм Ўтрассена
double* Str_alg(double* matrix1, double* matrix2, int N, int threshold)
{
	double* Rez;

	if (N <= threshold) // меньше критического значени€
		Rez = Simple_Mmult(matrix1, matrix2, N);
	else
	{
		Rez = CreateMatrix(N);
		N = N / 2;// разбили матрицу на 4 размеров Ќ/2 кажда€

		double* A[4]; double* B[4]; double* C[4]; double* P[7];

		double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
		double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

		// выделение пам€ти под вспомогательные матрицы 
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

		// выполнение рекурсивного умножени€ 
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

		// освобождение пам€ти 
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

double* Str_alg_tbb(double* matrix1, double* matrix2, int N, int threshold)
{
	double* Rez;

	if (N <= threshold)
		Rez = Simple_Mmult(matrix1, matrix2, N);
	else
	{
		Rez = CreateMatrix(N);
		N = N / 2;

		double* A[4]; double* B[4]; double* C[4]; double* P[7];

		double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
		double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

		for (int i = 0; i < 4; i++) {
			A[i] = CreateMatrix(N);
			B[i] = CreateMatrix(N);
		}


		parallel_for(blocked_range<size_t>(0, N),
			[=](const blocked_range<size_t>& r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					for (size_t j = 0; j < N; ++j)
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
				}
			});

		task_group TMP;
		TMP.run([&] { TMP1 = Add(A[0], A[3], N); });
		TMP.run([&] { TMP2 = Add(B[0], B[3], N); });
		TMP.run([&] { TMP3 = Add(A[2], A[3], N); });
		TMP.run([&] { TMP4 = Sub(B[1], B[3], N); });
		TMP.run([&] { TMP5 = Sub(B[2], B[0], N); });
		TMP.run([&] { TMP6 = Add(A[0], A[1], N); });
		TMP.run([&] { TMP7 = Sub(A[2], A[0], N); });
		TMP.run([&] { TMP8 = Add(B[0], B[1], N); });
		TMP.run([&] { TMP9 = Sub(A[1], A[3], N); });
		TMP.run([&] { TMP10 = Add(B[2], B[3], N); });
		TMP.wait();

		task_group p;
		p.run([&] { P[0] = Str_alg(TMP1, TMP2, N, threshold); });
		p.run([&] { P[1] = Str_alg(TMP3, B[0], N, threshold); });
		p.run([&] { P[2] = Str_alg(A[0], TMP4, N, threshold); });
		p.run([&] { P[3] = Str_alg(A[3], TMP5, N, threshold); });
		p.run([&] { P[4] = Str_alg(TMP6, B[3], N, threshold); });
		p.run([&] { P[5] = Str_alg(TMP7, TMP8, N, threshold); });
		p.run([&] { P[6] = Str_alg(TMP9, TMP10, N, threshold); });
		p.wait();

		task_group c;
		c.run([&] { C[0] = Sub(P[0], P[3], P[6], P[4], N); });
		c.run([&] { C[1] = Add(P[2], P[4], N); });
		c.run([&] { C[2] = Add(P[1], P[3], N); });
		c.run([&] { C[3] = Sub(P[0], P[2], P[5], P[1], N); });
		c.wait();

		parallel_for(blocked_range<size_t>(0, N),
			[=](const blocked_range<size_t>& r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					for (size_t j = 0; j < N; ++j)
					{
						Rez[i * 2 * N + j] = C[0][i * N + j];
						Rez[i * 2 * N + j + N] = C[1][i * N + j];
						Rez[i * 2 * N + j + 2 * N * N] = C[2][i * N + j];
						Rez[i * 2 * N + j + 2 * N * N + N] = C[3][i * N + j];
					}
				}
			});

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



int main()
{
	double* matr_A = NULL;
	double* matr_B = NULL;
	double* matr_Rez_Str = NULL;
	double* matr_Rez_Tbb = NULL;
	double StartStrAlg = 0;	// врем€ старта последовательного алгоритма Ўтрассена
	double TimeStrAlg = 0;// врем€, затраченное на работу последовательной версии алгоритма Ўтрассена
	int thr = 64;
	int k = 9;
	double StartTbb = 0;
	double TimeTbb = 0;

	int N = (int)pow(2.0, k);
	matr_A = CreateMatrix(N);
	matr_B = CreateMatrix(N);

	GenerateRandomMatrix(matr_A, matr_B, N);
	Null(matr_A, matr_B, N);

	// последовательный алгоритм 
	tick_count t1 = tick_count::now();

	matr_Rez_Str = Str_alg(matr_A, matr_B, N, thr);

	TimeStrAlg = (tick_count::now() - t1).seconds();

	cout << "Strassen linear:" << TimeStrAlg << endl;


	// tbb alg
	tick_count t0= tick_count::now();

	matr_Rez_Tbb = Str_alg_tbb(matr_A, matr_B, N, thr);
	TimeTbb = (tick_count::now() - t0).seconds();
	cout << "Strassen TBB: " << TimeTbb << endl;

	
	delete[] matr_Rez_Str;
	delete[] matr_Rez_Tbb;
	delete[] matr_B;
	delete[] matr_A;

	return 0;

}

