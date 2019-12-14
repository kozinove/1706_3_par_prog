// integral.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <math.h>
#include "mpi.h"

double eps = 0.00001;

// вычисляем двойной интеграл от этой функции 2х переменных:
double fun_xy(double x, double y) {
	return (cos(x) + cos(y));
}


// область интегрирования D
double fun_D(double x, double y) {
	double d;
	if ((x*x + y*y) < 1.0) d = 1.0;
	else d = 0.0;
	return d;
}


double calc_integral(double x1, double x2, double y1, double y2, int nx, int ny) {
	double dx = (x2 - x1) / nx;
	double dy = (y2 - y1) / ny;
	double sum = 0;
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {

			double x = x1 + ix * dx + dx * .5; // central
			double y = y1 + iy * dy + dy * .5; // central
			double f = fun_xy(x, y); // значение функции в точке
			double d = fun_D(x, y); // проверка, попадает ли точка в область интегрирования
			sum += f * d * dx * dy;
		}
	}
	return sum;
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int myrank, mysize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
	std::cout << " # rank: " << myrank << " size: " << mysize << std::endl;

	MPI_Status status;

	double x1, x2, y1, y2;
	double x1_loc, x2_loc, y1_loc, y2_loc;

	int nx, ny;
	if (myrank == 0) { // master process

		x1 = -1;
		x2 = 1;
		y1 = -1;
		y2 = 1;
		nx = 100;
		ny = 100;


		double dx = (x2 - x1) / (mysize - 1);

		for (int i = 1; i < mysize; i++) {
			x1_loc = x1 + dx * (i-1);
			x2_loc = x1_loc + dx;
			y1_loc = y1;
			y2_loc = y2;

			// send data to all slave processes (2-points exchange)
			MPI_Send(&x1_loc, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
			MPI_Send(&x2_loc, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
			MPI_Send(&y1_loc, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
			MPI_Send(&y2_loc, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);

			MPI_Send(&nx, 1, MPI_INTEGER, i, 99, MPI_COMM_WORLD);
			MPI_Send(&ny, 1, MPI_INTEGER, i, 99, MPI_COMM_WORLD);

		}
		double sum = 0, all_sum;
		MPI_Reduce(&sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		std::cout << " ------------------------- " << std::endl << std::endl;
		std::cout << " Integral = " << all_sum << std::endl;

		// check

		double check_sum = calc_integral(x1, x2, y1, y2, 1000, 1000);
		std::cout << " check = " << check_sum << std::endl;

	}

	// slave processes
	else {
		// get data on slave processes

		MPI_Recv(&x1_loc, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&x2_loc, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&y1_loc, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&y2_loc, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&nx, 1, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&ny, 1, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, &status);

		double sum1, sum2;
		sum1 = calc_integral(x1_loc, x2_loc, y1_loc, y2_loc, nx, ny);
		for (;;) {
			nx = nx * 2;
			ny = ny * 2;
			sum2 = calc_integral(x1_loc, x2_loc, y1_loc, y2_loc, nx, ny);
			if (abs(sum1 - sum2) > eps) {
				sum1 = sum2;
			}
			else break; // выход из бесконечного цикла
		}
		double all_sum = 0;
		MPI_Reduce(&sum2, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	// all processes


	MPI_Finalize();
	return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
