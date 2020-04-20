#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <time.h>
#include <tbb\parallel_reduce.h>
#include <tbb\blocked_range.h>
#include <tbb\task_scheduler_init.h>

#define Eps (0.000000000001)
#define MainProcess (0)
typedef time_t TimeWork;
int num_threads = 8;
double Serial_result = 0.;
double Parallel_result = 0.;

inline bool Test(double sequent, double parallel)
{
	bool OK = true;
	if (abs(sequent - parallel) > Eps)
		OK = false;

	std::cout << "\n Check the results ...";
	if (OK != true)
		std::cout << "Warning!!! Something went wrong." << std::endl;
	else
		std::cout << "Successfully!!!" << std::endl;

	return OK;
}

inline void CheckResults(double sequentTimeWork, double parallTimeWork)
{
	std::cout << "\n Who is faster? ...";
	if (parallTimeWork < sequentTimeWork)
		std::cout << " Parallel algorithm" << std::endl;
	else
		std::cout << " Sequential algorithm" << std::endl;

	std::cout.precision(3);
	std::cout.setf(std::ios::fixed);
	std::cout << " Speedup: " << sequentTimeWork / parallTimeWork << std::endl;
}

double Simpson_integration(double(*f)(double, double), double a, double b, int N, double fixX)
{
	double h = static_cast<double>(b - a) / N;

	double S4 = f(fixX, a + h), S2 = 0.;

	for (int i = 3; i < N; i += 2)
	{
		S4 += f(fixX, a + i * h);
		S2 += f(fixX, a + (i - 1) * h);
	}

	return (h / 3) * (f(fixX, a) + 4 * S4 + 2 * S2 + f(fixX, b));
}

class simpson_tbb_reduce
{
private:
	double a1, a2, b1, h2;
	int N;
	double(*F)(double, double);
	double I2;
	double I4;
public:
	explicit simpson_tbb_reduce(double(*_F)(double, double), double _a1, double _b1, int _N, double _a2, double _h2) : F(_F), a1(_a1), b1(_b1), N(_N), a2(_a2), h2(_h2), I2(0), I4(0) {}

	simpson_tbb_reduce(const simpson_tbb_reduce& m, tbb::split) : F(m.F), a1(m.a1), b1(m.b1), N(m.N), a2(m.a2), h2(m.h2), I2(0), I4(0) {}

	void operator()(const tbb::blocked_range<int>& r)
	{
		int end = r.end();
		for (int j = r.begin(); j < end; ++j)
		{
			if (j % 2 == 1)
				I4 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
			else
				I2 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
		}
	}

	void join(const simpson_tbb_reduce& simson_functor)
	{
		I2 += simson_functor.I2;
		I4 += simson_functor.I4;
	}

	double getI2() { return I2; }
	double getI4() { return I4; }
};


TimeWork Simpson2(double(*F)(double, double), double a1, double b1, double a2, double b2, int n)
{
	time_t time_start, time_end;

	time_start = clock();  //time(&time_start);

	if (n <= 0)
		return 0.;

	int sign = 1;

	if (abs(a1 - b1) < Eps || abs(a2 - b2) < Eps)
		return 0.;

	if (a1 > b1)
	{
		double tmp = b1;
		b1 = a1;
		a1 = b1;
		sign *= -1;
	}

	if (a2 > b2)
	{
		double tmp = b2;
		b2 = a2;
		a2 = b2;
		sign *= -1;
	}

	int N = 2 * n;

	double h2 = static_cast<double>(b2 - a2) / N;

	double I2 = 0., I4 = 0.;

	for (int j = 1; j < N; ++j)
	{
		if (j % 2 == 1)
			I4 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
		else
			I2 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
	}

	Serial_result = sign * (h2 / 3) * (Simpson_integration(F, a1, b1, N, a2) + 2 * I2 + 4 * I4 + Simpson_integration(F, a1, b1, N, b2));

	time_end = clock();  //time(&time_end);
	return (time_end - time_start);
}

double Simpson2_OMP_parallel(double(*F)(double, double), double a1, double b1, double a2, double b2, int n)
{
	if (n <= 0)
		return 0.;

	int sign = 1;

	if (abs(a1 - b1) < Eps || abs(a2 - b2) < Eps)
		return 0.;

	if (a1 > b1)
	{
		double tmp = b1;
		b1 = a1;
		a1 = b1;
		sign *= -1;
	}

	if (a2 > b2)
	{
		double tmp = b2;
		b2 = a2;
		a2 = b2;
		sign *= -1;
	}

	int N = 2 * n;

	double h2 = static_cast<double>(b2 - a2) / N;

	double I2 = 0., I3 = 0.;

	int j = 0;

#pragma omp parallel for shared(a1, b1, a2, N, h2) private(j) reduction(+:I3,I2)
	
	for (j = 1; j < N; ++j)
	{
		if (j % 2 == 1)
			I3 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
		else
			I2 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
	}

	return sign * (h2 / 3) * (Simpson_integration(F, a1, b1, N, a2) + 2 * I2 + 4 * I3 + Simpson_integration(F, a1, b1, N, b2));
}

double Simpson2_TBB_parallel(double(*F)(double, double), double a1, double b1, double a2, double b2, int n)
{
	if (n <= 0)
		return 0.;

	int sign = 1;

	if (abs(a1 - b1) < Eps || abs(a2 - b2) < Eps)
		return 0.;

	if (a1 > b1)
	{
		double tmp = b1;
		b1 = a1;
		a1 = b1;
		sign *= -1;
	}

	if (a2 > b2)
	{
		double tmp = b2;
		b2 = a2;
		a2 = b2;
		sign *= -1;
	}

	int N = 2 * n;

	double h2 = static_cast<double>(b2 - a2) / N;

	tbb::task_scheduler_init init(4);
	simpson_tbb_reduce simpson_tbb_body(F, a1, b1, N, a2, h2);

	tbb::parallel_reduce(tbb::blocked_range<int>(1, N, 5), simpson_tbb_body);
	init.terminate();

	return sign * (h2 / 3) * (Simpson_integration(F, a1, b1, N, a2) + 2 * simpson_tbb_body.getI2() + 4 * simpson_tbb_body.getI4() + Simpson_integration(F, a1, b1, N, b2));
}


inline double func(double x, double y)
{
	return (x * x * y - x) / sqrt(x + y * y * x);
}


inline void Initialize_from_Command_line(int argc, char* argv[], int* N, int* num_threads)
{
	if (argc > 1)
	{
		*num_threads = atoi(argv[1]);
		*N = (argc > 2) ? atoi(argv[2]) : 1000;
	}

	if (*N <= 0 || *N > 1000000000)
		*N = 1000;

	if (*num_threads <= 0 || *num_threads > 50)
		*num_threads = 2;
}



int main(int argc, char* argv[])
{
	TimeWork parallTimeWork = 0., sequentTimeWork = 0.;
	time_t time_start, time_end;
	int ProcNum = 0;
	int	ProcRank = 0;
	int n = 10000;


	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	//MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	Initialize_from_Command_line(argc, argv, &n, &num_threads);

	//if (ProcRank == MainProcess)
	//{
	    time(&time_start);
		sequentTimeWork = Simpson2(&func, 1, 2, 2, 3, n);
		time(&time_end);
		std::cout.setf(std::ios::fixed);
		std::cout.precision(15);
		std::cout << "\n ******  Simpson Multi-Dimensional Integration ******";
		std::cout << "\n ****** " << "num_processes/num_threads = " << ProcNum << " / " << num_threads << " ******";
		std::cout << "\n ****** " << "num_iterations = " << n << " ******\n";

		std::cout << "\nSimpson integrate = " << Serial_result << std::endl;
		std::cout.precision(6);
		std::cout << " time = " << sequentTimeWork << "ms\n ***************************";
	//}
	time_start = clock();  //time(&time_start);
	Parallel_result = Simpson2_TBB_parallel(&func, 1, 2, 2, 3, n);
	time_end = clock();  //time(&time_end);
	parallTimeWork = (time_end - time_start);
	//if (ProcRank == MainProcess)
	//{
	std::cout << "\n ******  TBB  ******";
		std::cout.precision(15);
		std::cout << "\nSimpson integrate = " << Parallel_result;
		std::cout.precision(6);
		if (Test(Serial_result, Parallel_result))
			std::cout << " time = " << parallTimeWork << "ms\n ***************************";

		CheckResults(sequentTimeWork, parallTimeWork);
		std::cout << "\n ***************************\n";
	//}

		time_start = clock();//time(&time_start); 
		Parallel_result = Simpson2_OMP_parallel(&func, 1, 2, 2, 3, n);
		time_end = clock(); // time(&time_end);
		parallTimeWork = (time_end - time_start);
		//if (ProcRank == MainProcess)
		//{
		std::cout << "\n ******  OpenMP  ******";
		std::cout.precision(15);
		std::cout << "\nSimpson integrate = " << Parallel_result;
		std::cout.precision(6);
		if (Test(Serial_result, Parallel_result))
			std::cout << " time = " << parallTimeWork << "ms\n ***************************";

		CheckResults(sequentTimeWork, parallTimeWork);
		std::cout << "\n ***************************\n";
		//}

	//MPI_Finalize();
	return 0;
}
