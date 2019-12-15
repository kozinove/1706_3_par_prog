#include <mpi.h>
#include <iostream>


#define ESP (0.000000000001)
#define ROOTPROCESS (0)
typedef double TimeWork;
double Serial_result = 0.;
double Parallel_result = 0.;

inline bool Test(double sequent, double parallel)
{
	bool OK = true;
	if (abs(sequent - parallel) > ESP)
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

TimeWork Simpson2(double(*F)(double, double), double a1, double b1, double a2, double b2, int n)
{
	double time_start, time_end;

	time_start = MPI_Wtime();

	if (n <= 0)
		return 0.;

	int sign = 1;

	if (abs(a1 - b1) < ESP || abs(a2 - b2) < ESP)
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

	for (int j = 1; j < N; ++j)
	{
		if (j % 2 == 1)
			I3 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
		else
			I2 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
	}

	Serial_result = sign * (h2 / 3) * (Simpson_integration(F, a1, b1, N, a2) + 2 * I2 + 4 * I3 + Simpson_integration(F, a1, b1, N, b2));

	time_end = MPI_Wtime();
	return (time_end - time_start) * 1000;
}

TimeWork Simpson2_MPI_parallel(double(*F)(double, double), double a1, double b1, double a2, double b2, int n, int ProcRank, int ProcNum)
{
	double time_start, time_end;

	MPI_Barrier(MPI_COMM_WORLD);
	time_start = MPI_Wtime();

	if (n <= 0)
		return 0.;

	int sign = 1;

	if (abs(a1 - b1) < ESP || abs(a2 - b2) < ESP)
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


	int k = (N - 1) / ProcNum;
	int m = (N - 1) % ProcNum;

	int diapason0 = 1;
	int diapason1 = k + 1 + ((m != 0) ? 1 : 0);

	for (int i = 1; i < (ProcRank + 1); ++i)
	{
		diapason0 = diapason1;
		diapason1 += k + ((i < m) ? 1 : 0);
	}


	double local_I2 = 0.;
	double local_I3 = 0.;
	int j = 0;

	for (j = diapason0; j < diapason1; ++j)
	{
		if (j % 2 == 1)
			local_I3 += Simpson_integration(F, a1, b1, N, a2 + j * h2);
		else
			local_I2 += Simpson_integration(F, a1, b1, N, a2 + j * h2);

	}

	double I2 = 0., I3 = 0.;

	MPI_Reduce(&local_I2, &I2, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&local_I3, &I3, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);


	Parallel_result = sign * (h2 / 3) * (Simpson_integration(F, a1, b1, N, a2) + 2 * I2 + 4 * I3 + Simpson_integration(F, a1, b1, N, b2));

	time_end = MPI_Wtime();
	return (time_end - time_start) * 1000;
}

inline double func2(double x, double y)
{
	return x * y * y;
}

inline double func1(double x, double y)
{
	return (x * x * y - x) / sqrt(x + y * y * x);
}


inline void Initialize_from_Command_line(int argc, char* argv[], int* N)
{
	if (argc > 1)
	{
		*N = (argc > 1) ? atoi(argv[1]) : 1000;
	}

	if (*N <= 0 || *N > 1000000000)
		*N = 1000;
}



int main(int argc, char* argv[])
{
	TimeWork parallTimeWork = 0., sequentTimeWork = 0.;
	int ProcNum; // Количество процессов
	int	ProcRank; // Ранг процесса
	int n = 1000;
	Initialize_from_Command_line(argc, argv, &n);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == ROOTPROCESS)
	{
		sequentTimeWork = Simpson2(&func1, 1, 2, 2, 3, n);

		std::cout.setf(std::ios::fixed);
		std::cout.precision(15); 
		std::cout << "\n ******  Simpson Integration  ******";
		std::cout << "\nnum_processes = " << ProcNum;
		std::cout << "\nnum_iterations = " << 2 * n << "\n";

		std::cout << "\nSimpson integrate = " << Serial_result << std::endl;
		std::cout.precision(6);
		std::cout << " time = " << sequentTimeWork << "ms\n";
	}

	parallTimeWork = Simpson2_MPI_parallel(&func1, 1, 2, 2, 3, n, ProcRank, ProcNum);

	if (ProcRank == ROOTPROCESS)
	{
		std::cout.precision(15);
		std::cout << "\nSimpson integrate = " << Parallel_result;
		std::cout.precision(6);
		if (Test(Serial_result, Parallel_result))
			std::cout << " time = " << parallTimeWork << "ms\n";

		CheckResults(sequentTimeWork, parallTimeWork);
		std::cout << "\n ***************************\n";
	}

	return MPI_Finalize();
}
