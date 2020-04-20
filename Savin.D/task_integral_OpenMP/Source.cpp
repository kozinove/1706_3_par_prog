#include <iostream>
#include <math.h>
#include <omp.h>

double eps = 0.00001;

// вычисляем двойной интеграл от этой функции 2х переменных:
double fun_xy(double x, double y) {
	return (cos(x) + cos(y));
}


// область интегрирования D
double fun_D(double x, double y) {
	double d;
	if ((x * x + y * y) < 1.0) d = 1.0;
	else d = 0.0;
	return d;
}



double calc_integral_1D_openmp(double x, double y1, double y2, int ny)
{

	double sum = 0;
	double dy = (y2 - y1) / ny;
	int iy = 0;

#pragma omp parallel reduction(+:sum)
	{
#pragma omp for schedule(static)
		for (iy = 0; iy < ny; iy++)
		{
			double y = y1 + iy * dy + dy * .5; // central
			double f = fun_xy(x, y); // значение функции в точке
			double d = fun_D(x, y); // проверка, попадает ли точка в область интегрирования
			sum += f * d * dy;
		}
	}
	return sum;
}

double calc_integral_openmp(double x1, double x2, double y1, double y2, int nx, int ny) {
	double dx = (x2 - x1) / nx;
	double sum = 0;
	for (int ix = 0; ix < nx; ix++) {


		double x = x1 + ix * dx + dx * .5; // central
		double integral_1D = calc_integral_1D_openmp(x, y1, y2, ny);
		sum += integral_1D * dx;

	}
	return sum;
}

double calc_integral_1D(double x, double y1, double y2, int ny) {
	double sum = 0;
	double dy = (y2 - y1) / ny;
	for (int iy = 0; iy < ny; iy++) {
		double y = y1 + iy * dy + dy * .5; // central
		double f = fun_xy(x, y); // значение функции в точке
		double d = fun_D(x, y); // проверка, попадает ли точка в область интегрирования
		sum += f * d * dy;
	}
	return sum;
}


double calc_integral(double x1, double x2, double y1, double y2, int nx, int ny) {
	double dx = (x2 - x1) / nx;
	double sum = 0;
	for (int ix = 0; ix < nx; ix++) {


		double x = x1 + ix * dx + dx * .5; // central
		double integral_1D = calc_integral_1D(x, y1, y2, ny);
		sum += integral_1D * dx;

	}
	return sum;
}


int main(int argc, char** argv) {


	double x1, x2, y1, y2;
	x1 = -1;
	x2 = 1;
	y1 = -1;
	y2 = 1;

	int thread;
	std::cout << "thread = ";
	std::cin >> thread;
	omp_set_num_threads(thread);

	double start = 0;
	double end = 0;
	double start_openmp = 0;
	double end_openmp = 0;
     

	start = omp_get_wtime();
	double check_sum = calc_integral(x1, x2, y1, y2, 1000, 1000);
	end = omp_get_wtime();
	std::cout << " posledovatelno check = " << check_sum << std::endl;

	start_openmp = omp_get_wtime();
	double check_sum_openmp = calc_integral_openmp(x1, x2, y1, y2, 1000, 1000);
	end_openmp = omp_get_wtime();
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "paralel check = " << check_sum_openmp << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "posledovatelno = " << end - start << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "paralel = " << end_openmp - start_openmp << std::endl;

	return 0;

}


