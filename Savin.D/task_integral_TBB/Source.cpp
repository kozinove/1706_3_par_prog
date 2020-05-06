#include <iostream>
#include <math.h>
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"
#include "tbb/parallel_reduce.h"
#include "tbb/tbb.h"

using namespace std;
using namespace tbb;

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


//
//double calc_integral_1D_openmp(double x, double y1, double y2, int ny)
//{
//
//	double sum = 0;
//	double dy = (y2 - y1) / ny;
//	int iy = 0;

//#pragma omp parallel reduction(+:sum)
//	{
//#pragma omp for schedule(static)
//		for (iy = 0; iy < ny; iy++)
//		{
//			double y = y1 + iy * dy + dy * .5; // central
//			double f = fun_xy(x, y); // значение функции в точке
//			double d = fun_D(x, y); // проверка, попадает ли точка в область интегрирования
//			sum += f * d * dy;
//		}
//	}
//	return sum;
//}
//
//double calc_integral_openmp(double x1, double x2, double y1, double y2, int nx, int ny) {
//	double dx = (x2 - x1) / nx;
//	double sum = 0;
//	for (int ix = 0; ix < nx; ix++) {
//
//
//		double x = x1 + ix * dx + dx * .5; // central
//		double integral_1D = calc_integral_1D_openmp(x, y1, y2, ny);
//		sum += integral_1D * dx;
//
//	}
//	return sum;
//}

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


double calc_integral_tbb(double x1, double x2, double y1, double y2, int nx, int ny) {
	double dx = (x2 - x1) / nx;
	double sum = 0;
	parallel_for (blocked_range<int>(0, nx),
		[&](const blocked_range<int>& r)
		{
			for (int ix = r.begin(); ix != r.end(); ++ix)
			{
				double x = x1 + ix * dx + dx * .5; // central
				double integral_1D = calc_integral_1D(x, y1, y2, ny);
				sum += integral_1D * dx;
			}
		}
	);
	return sum;
}

double calc_integral_tbb1(double x1, double x2, double y1, double y2, int nx, int ny) {
	double dx = (x2 - x1) / nx;
	double sum = 0;
	return (parallel_for(blocked_range<int>(0, nx),
		[&](const blocked_range<int>& r, double sum1 = 0)
		{
			for (int ix = r.begin(); ix != r.end(); ++ix)
			{
				double x = x1 + ix * dx + dx * .5; // central
				double integral_1D = calc_integral_1D(x, y1, y2, ny);
				sum1 += integral_1D * dx;
			}
			return sum1;
		}
	));
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
	task_scheduler_init init(thread);

	double time = 0;
	double time_tbb = 0;


	tick_count t0 = tick_count::now();
	double check_sum = calc_integral(x1, x2, y1, y2, 1000, 1000);
	time = (tick_count::now() - t0).seconds();
	std::cout << " posledovatelno check = " << check_sum << std::endl;

	tick_count t1 = tick_count::now();
	double check_sum_tbb = calc_integral_tbb(x1, x2, y1, y2, 1000, 1000);
	time_tbb = (tick_count::now() - t1).seconds();
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "paralel check = " << check_sum_tbb << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "posledovatelno = " << time << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "paralel = " << time_tbb << std::endl;

	return 0;

}


