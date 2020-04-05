#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>

using namespace std;

void Gauss(double* A, double* b, double* x, int size) 
{
	int i, j, k;
	double alfa;
	double* matrix = new double[size * size];
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			matrix[i * size + j] = A[i * size + j];
	double* vector = new double[size];
	for (int i = 0; i < size; i++)
		vector[i] = b[i];
	
	for (j = 0; j < size; j++)
	{
		for (i = j + 1; i < size; i++)
		{
			alfa = matrix[i*size + j] / matrix[j*size + j];
			for (k = j; k < size; k++)
			{
				matrix[i * size + k] -= alfa * matrix[j * size + k];
			}
			vector[i] -= alfa * vector[j];
		}
	}

	x[size - 1] = vector[size - 1] / matrix[(size - 1) * size + size - 1];
	for (i = size - 2; i >= 0; i--)
	{
		double sum = 0;
		for (j = i + 1; j < size; j++)
		{
			sum += matrix[i * size + j] * x[j];
		}
		x[i] = (vector[i] - sum) / matrix[i * size + i];
	}
	/*cout << "Gauss result: " << endl;
	for (i = 0; i < size; i++)
	{
		cout << x[i] << endl;
	}*/

	delete[] matrix;
	delete[] vector;
}
void GenerateStartData(double* &Matrix, double*& rightVector ,const int &size, int argc, char* argv[],const int &ProcRank,const int &ProcNum, MPI_Comm COMM_WORKING)
{
	int max = 0;
	int min = 0;
	Matrix = new double[size * size];
	rightVector = new double[size];
	int* sendCount = new int[ProcNum];
	int* displs = new int[ProcNum];
	int startNum = size / ProcNum;
	double* transpMatrix = nullptr;
	double* transpBuffer = nullptr;
	for (int i = 0; i < ProcNum; i++)
	{
		sendCount[i] = startNum * size;
		displs[i] = 0;
	}
	for (int i = 0; i < size - startNum * ProcNum; i++)
		sendCount[i] += size;
	for (int i = 1; i < ProcNum; i++)
		displs[i] = displs[i - 1] + sendCount[i - 1];

	if (ProcRank == 0)
	{
		min = 0;
		max = 10;
		if (argc > 2)
			min = atoi(argv[2]);
		if (argc > 3)
			max = atoi(argv[3]);
		if (max == min)
		{
			if (ProcRank == 0)
				cout << "ERROR DATA, WILL BE USED MAX = 10, MIN = 0" << endl;
			max = 10;
			min = 0;
		}
		else if (max < min)
		{
			if (ProcRank == 0)
				cout << "ERROR DATA, MIN AND MAX WILL BE SWAPPED" << endl;
			swap(min, max);
		}
		int dopIter = 0;
		if (argc > 4)
		{
			dopIter = atoi(argv[4]);
		}
		cout << "USING DATA: " << endl << "Size = " << size << endl << "Approximate min = " << min << endl << "Approximate max = " << max << endl << "Additional iterations = " << dopIter << endl;;

		for (int i = 0; i < size * size; i++)
			Matrix[i] = 0;
		srand(time(NULL));
		transpMatrix = new double[size * size];
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				Matrix[i * size + j] = (double)(rand()) / RAND_MAX * (max - min) + min;
				//Matrix[i * size + j] = rand() % 11;
			}
			rightVector[i] = ((double)(rand()) / RAND_MAX * (max - min) + min);
			//rightVector[i] = rand() % 11;
		}

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				transpMatrix[i * size + j] = Matrix[j * size + i];
			}


	}
	MPI_Bcast(Matrix, size * size, MPI_DOUBLE, 0, COMM_WORKING);
	MPI_Bcast(rightVector, size, MPI_DOUBLE, 0, COMM_WORKING);
	MPI_Bcast(&max, 1, MPI_INT, 0, COMM_WORKING);
	transpBuffer = new double[sendCount[ProcRank]];
	for (int i = 0; i < size; i++)
		transpBuffer[i] = 0;
	MPI_Scatterv(transpMatrix, sendCount, displs, MPI_DOUBLE, transpBuffer, sendCount[ProcRank], MPI_DOUBLE, 0, COMM_WORKING);




	double* resBuffer = new double[sendCount[ProcRank]];
	for (int i = 0; i < sendCount[ProcRank]; i++)
		resBuffer[i] = 0;


	for (int i = 0; i < sendCount[ProcRank] / size; i++)
		for (int j = 0; j < size; j++)
			for (int k = 0; k < size; k++)
				resBuffer[i * size + j] += transpBuffer[i * size + k] * Matrix[j + k * size];

	for (int i = 0; i < sendCount[ProcRank]; i++)
		resBuffer[i] = resBuffer[i] / (max * size);
	
	MPI_Barrier(COMM_WORKING);
	MPI_Gatherv(resBuffer, sendCount[ProcRank], MPI_DOUBLE, Matrix, sendCount, displs, MPI_DOUBLE, 0, COMM_WORKING);

	if (ProcRank == 0)
	{
		/*cout << "----------------------" << endl;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				cout << Matrix[i * size + j] << "  ";
			cout << endl;
		}*/
		delete[] transpMatrix;
	}
	else
		delete[] Matrix;
	delete[] transpBuffer;
	delete[] resBuffer;
	delete[] sendCount;
	delete[] displs;

}
double* MmultipleV(double* matrix, double* vector, const int& size) {
	double* result = new double[size];
	for (int i = 0; i < size; i++) {
		result[i] = 0;
		for (int j = 0; j < size; j++)
			result[i] += matrix[size * i + j] * vector[j];
	}
	return result;
}
double* VminusV(double* vectorLeft, double* vectorRight, const int& size) {
	double* result = new double[size];
	for (int i = 0; i < size; i++) {
		result[i] = vectorLeft[i] - vectorRight[i];
	}
	return result;
}
double* VplusV(double* vectorLeft, double* vectorRight, const int& size) {
	double* result = new double[size];
	for (int i = 0; i < size; ++i) {
		result[i] = vectorLeft[i] + vectorRight[i];
	}
	return result;
}
double ScalMul(double* vectorLeft, double* vectorRight, const int& size) {
	double result = 0;
	for (int i = 0; i < size; ++i) {
		result += vectorLeft[i] * vectorRight[i];
	}
	return result;
}
double* VmultipleNum(double* vector, double num, const int& size) {
	double* result = new double[size];
	for (int i = 0; i < size; ++i) {
		result[i] = vector[i] * num;
	}
	return result;
}

void ConsistentGradient(double* matrix, double* rightVector, double* &result, int size)
{
	int iter = 0;
	
	double* prevResult = new double[size];
	for (int i = 0; i < size; i++)
		prevResult[i] = 0;
	double* r = new double[size];
	double* prevR = new double[size];
	double* p = new double[size];

	r = VminusV(rightVector, MmultipleV(matrix, prevResult, size), size);
	//cout << "r  ";
	for (int i = 0; i < size; i++)
	{
		//cout << r[i] << "  ";
		p[i] = r[i];
	}
	cout << endl;
	double alfa = 0;
	double beta = 0;
	while (iter < size)
	{
		/*cout << "Iter = " << iter << endl;
		for (int i = 0; i < size; i++)
		{
			cout << "prevResult: " << prevResult[i] << endl;
		}*/
		alfa = ScalMul(r, r, size) / ScalMul(MmultipleV(matrix, p, size), p, size);
		//cout << "alfa" << alfa << endl;
		result = VplusV(prevResult, VmultipleNum(p, alfa, size), size);
		/*for (int i = 0; i < size; i++)
		{
			cout << "p" << p[i] << endl;
		}*/
		/*for (int i = 0; i < size; i++)
		{
			cout << "poslRes" << result[i] << endl;
		}*/
		for (int i = 0; i < size; i++)
		{
			prevR[i] = r[i];
		}
		r = VminusV(r, VmultipleNum(MmultipleV(matrix, p, size), alfa, size), size);
		/*for (int i = 0; i < size; i++)
		{
			cout << "r" << r[i] << endl;
		}*/
		beta = ScalMul(r, r, size) / ScalMul(prevR, prevR, size);
		//cout << "sm1:" << ScalMul(r, r, size) << " sm2:" << ScalMul(prevR, prevR, size) << " beta: " << beta << endl;
		//cout << " beta: " << beta << endl;
		p = VplusV(r, VmultipleNum(p, beta, size), size);
		/*for (int i = 0; i < size; i++)
		{
			cout << "p" << p[i] << endl;
		}*/
		iter++;

		for (int i = 0; i < size; i++)
			prevResult[i] = result[i];
	}
	/*cout << "Consistent gradient result: " << endl;
	for (int i = 0; i < size; i++)
		cout << result[i] << endl;*/
}

void main(int argc, char* argv[])
{
	int ProcNum = 0;
	int ProcRank = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	if (argc < 5)
	{
		if (ProcRank == 0)
			cout << "NO ENOUGH ARGUMENTS, WILL BE USED SERVICE DATA" << endl;
	}
	if (argc > 6)
	{
		if (ProcRank == 0)
			cout << "TOO MUCH ARGUMENTS, SOME WON`T BE USED" << endl;
	}

	double* matrix = nullptr;
	double* rightVector = nullptr;
	int size = 5;
	if (argc > 1)
	{
		size = atoi(argv[1]);
	}
	if (size < 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED SIZE = |SIZE|" << endl;
		size = abs(size);
	}
	else if (size == 0)
	{
		if (ProcRank == 0)
			cout << "ERROR DATA, WILL BE USED SIZE = 5" << endl;
		size = 5;
	}

	MPI_Group Working;
	MPI_Comm_group(MPI_COMM_WORLD, &Working);

	int workingP;

	if (ProcNum > size)
		workingP = size;
	else workingP = ProcNum;
	int* arrWorkingP = new int[workingP];
	for (int i = 0; i < workingP; i++)
	{
		arrWorkingP[i] = i;
	}
	MPI_Group_incl(Working, workingP, arrWorkingP, &Working);
	delete[] arrWorkingP;
	MPI_Comm COMM_WORKING;
	MPI_Comm_create(MPI_COMM_WORLD, Working, &COMM_WORKING);
	MPI_Group_free(&Working);
	ProcNum = workingP;

	if (ProcRank < ProcNum)
	{
		GenerateStartData(matrix, rightVector, size, argc, argv, ProcRank, ProcNum, COMM_WORKING);
	}

	double* resultG = nullptr;
	double* result = new double[size];
	double startTime[3] = { 0 };
	double finishTime[3] = { 0 };
	double workTime[3] = { 0 };

	if (ProcRank == 0)
	{

		cout << "----------------------" << endl;
		resultG = new double[size];
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				cout << matrix[i * size + j] << "  ";
			cout << "       " << rightVector[i] << endl;
		}
		startTime[0] = MPI_Wtime();
		Gauss(matrix, rightVector, resultG, size);
		finishTime[0] = MPI_Wtime();
		workTime[0] = finishTime[0] - startTime[0];

		startTime[1] = MPI_Wtime();
		ConsistentGradient(matrix, rightVector, result, size);
		finishTime[1] = MPI_Wtime();
		workTime[1] = finishTime[1] - startTime[1];
		/*int checkCorr = 0;
		for (int i = 0; i < size; i++)
			checkCorr += result[i] - resultG[i];
		cout << "checkCorr = " << checkCorr << endl;
		if (checkCorr == 0)
			cout << "Consistent gradient result is correct!" << endl;
		else
			cout << "Wrong consistent gradient result!" << endl;*/
	}

	

	if (ProcRank < ProcNum)
	{
		int iter = 0;

		double* r = new double[size];
		double* prevR = new double[size];
		double* p = new double[size];
		double* prevResult = new double[size];
		for (int i = 0; i < size; i++)
			prevResult[i] = 0;
		for (int i = 0; i < size; i++)
			result[i] = 0;
		int* sendCount = new int[ProcNum];
		int* displs = new int[ProcNum];
		int* recvCount = new int[ProcNum];
		int* gatherDispls = new int[ProcNum];
		int startNum = size / ProcNum;
		for (int i = 0; i < ProcNum; i++)
		{
			sendCount[i] = startNum * size;
			displs[i] = 0;
		}
		for (int i = 0; i < size - startNum * ProcNum; i++)
			sendCount[i] += size;
		for (int i = 1; i < ProcNum; i++)
			displs[i] = displs[i - 1] + sendCount[i - 1];

		for (int i = 0; i < ProcNum; i++)
		{
			recvCount[i] = sendCount[i] / size;
			gatherDispls[i] = displs[i] / size;
		}

		/*if (ProcRank == 0)
		{
			cout << "SendCount: ";
			for (int i = 0; i < ProcNum; i++)
				cout << sendCount[i] << "  ";
			cout << endl << "Displs: ";
			for (int i = 0; i < ProcNum; i++)
				cout << displs[i] << "  ";
			cout << endl << "RecvCount: ";
			for (int i = 0; i < ProcNum; i++)
				cout << recvCount[i] << "  ";
			cout << endl << "GatherDispls: ";
			for (int i = 0; i < ProcNum; i++)
				cout << gatherDispls[i] << "  ";
			cout << endl;
		}*/
		double alfa = 0;
		double beta = 0;
		double scalMulR = 0;
		double scalMulMP = 0;
		double tempScal = 0;

		double* matrixBuffer = new double[sendCount[ProcRank]];
		double* resBuffer = new double[recvCount[ProcRank]];
		double* tempResult = new double[recvCount[ProcRank]];
		startTime[2] = MPI_Wtime();
		
		for (int i = 0; i < size; i++)
		{
			r[i] = rightVector[i];
			p[i] = r[i];
		}
		int dopIter = 0;
		if (argc > 4)
		{
			dopIter = atoi(argv[4]);
		}
		while (iter < size + dopIter)
		{
			tempScal = 0;
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempScal += r[i + gatherDispls[ProcRank]] * r[i + gatherDispls[ProcRank]];
			}
			//MPI_Barrier(COMM_WORKING);
			MPI_Allreduce(&tempScal, &scalMulR, 1, MPI_DOUBLE, MPI_SUM, COMM_WORKING);


			/*MPI_Barrier(COMM_WORKING);
			if (ProcRank == 0)
			{
				cout << "p Iter = " << iter << endl;
				for (int i = 0; i < size; i++)
				{
					cout << "p prevResult: " << prevResult[i] << endl;
				}
			}*/


			MPI_Scatterv(matrix, sendCount, displs, MPI_DOUBLE, matrixBuffer, sendCount[ProcRank], MPI_DOUBLE, 0, COMM_WORKING);

			for (int i = 0; i < recvCount[ProcRank]; i++)
				resBuffer[i] = 0;
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				for (int k = 0; k < size; k++)
					resBuffer[i] += matrixBuffer[i * size + k] * p[k];
			}
			tempScal = 0;
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempScal += resBuffer[i] * p[i + gatherDispls[ProcRank]];

			}
			//MPI_Barrier(COMM_WORKING);
			MPI_Allreduce(&tempScal, &scalMulMP, 1, MPI_DOUBLE, MPI_SUM, COMM_WORKING);

			alfa = scalMulR / scalMulMP;
			/*MPI_Barrier(COMM_WORKING);
			if (ProcRank == 0)
			{
				cout << "p alfa" << alfa << endl;
			}*/


			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempResult[i] = prevResult[i + gatherDispls[ProcRank]] + alfa * p[i + gatherDispls[ProcRank]];
			}
			MPI_Allgatherv(tempResult, recvCount[ProcRank], MPI_DOUBLE, result, recvCount, gatherDispls, MPI_DOUBLE, COMM_WORKING);

			/*MPI_Barrier(COMM_WORKING);
			if (ProcRank == 0)
			{
				for (int i = 0; i < size; i++)
				{
					cout << "p poslRes" << result[i] << endl;
				}
			}*/
			for (int i = 0; i < size; i++)
			{
				prevR[i] = r[i];
			}
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempResult[i] = r[i + gatherDispls[ProcRank]] - alfa * resBuffer[i];
			}
			MPI_Allgatherv(tempResult, recvCount[ProcRank], MPI_DOUBLE, r, recvCount, gatherDispls, MPI_DOUBLE, COMM_WORKING);
			/*MPI_Barrier(COMM_WORKING);
			if (ProcRank == 0)
			{
				for (int i = 0; i < size; i++)
				{
					cout << "p r" << r[i] << endl;
				}
			}*/
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempResult[i] = r[i + gatherDispls[ProcRank]] - alfa * resBuffer[i];
			}
			tempScal = 0;
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempScal += r[i + gatherDispls[ProcRank]] * r[i + gatherDispls[ProcRank]];
			}
			//MPI_Barrier(COMM_WORKING);
			MPI_Allreduce(&tempScal, &scalMulR, 1, MPI_DOUBLE, MPI_SUM, COMM_WORKING);
			tempScal = 0;
			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempScal += prevR[i + gatherDispls[ProcRank]] * prevR[i + gatherDispls[ProcRank]];
			}
			//MPI_Barrier(COMM_WORKING);
			MPI_Allreduce(&tempScal, &scalMulMP, 1, MPI_DOUBLE, MPI_SUM, COMM_WORKING);
			beta = scalMulR / scalMulMP;

			/*MPI_Barrier(COMM_WORKING);
			if (ProcRank == 0)
			{
					cout << "p beta" << beta << endl;
			}*/

			for (int i = 0; i < recvCount[ProcRank]; i++)
			{
				tempResult[i] = r[i + gatherDispls[ProcRank]] + beta * p[i + gatherDispls[ProcRank]];
			}
			MPI_Allgatherv(tempResult, recvCount[ProcRank], MPI_DOUBLE, p, recvCount, gatherDispls, MPI_DOUBLE, COMM_WORKING);


			for (int i = 0; i < size; i++)
				prevResult[i] = result[i];


			iter++;
			//MPI_Barrier(COMM_WORKING);

		}
		MPI_Barrier(COMM_WORKING);
		finishTime[2] = MPI_Wtime();
		workTime[2] = finishTime[2] - startTime[2];


		if (ProcRank == 0)
		{
			int checkCorr = 0;
			for (int i = 0; i < size; i++)
				checkCorr += abs(result[i] - resultG[i]);
			
			cout << "checkCorr = " << checkCorr << endl;
			if (checkCorr != 0)
			{
				cout << "Parallel gradient result       Gauss result          Error" << endl;
				for (int i = 0; i < size; i++)
					cout << result[i] << "                       " << resultG[i] << "                   " << result[i] - resultG[i] << endl;
			}
			else
			{
				cout << "Parallel gradient result: " << endl;
				for (int i = 0; i < size; i++)
					cout << result[i] << endl;
			}
			if (checkCorr == 0)
				cout << "Parallel gradient result is full correct!" << endl;
			else
				cout << "Parallel gradient result counted with error!" << endl;


			cout << "Time: " << endl;
			cout << "Gauss: " << workTime[0] << " sec." << endl;
			cout << "Consistent gradient: " << workTime[1] << " sec." << endl;
			cout << "Parallel gradient: " << workTime[2] << " sec." << endl;


		}
		MPI_Comm_free(&COMM_WORKING);
	}

	
	MPI_Finalize();
}


/*double* matrixBuffer = new double[sendCount[ProcRank]];
MPI_Scatterv(matrix, sendCount, displs, MPI_DOUBLE, matrixBuffer, sendCount[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

double* resBuffer = new double[sendCount[ProcRank]];

for (int i = 0; i < sendCount[ProcRank] / size; i++)
	for (int j = 0; j < size; j++)
		for (int k = 0; k < size; k++)
			resBuffer[i * size + j] += matrixBuffer[i * size + k] * prevResult[k];

for (int i = 0; i < sendCount[ProcRank]; i++)
{
	cout << ProcRank << ":  " << rightVector[displs[ProcRank] / size + i] << endl;
	resBuffer[i] = rightVector[displs[ProcRank] / size + i] - resBuffer[i];
	cout << ProcRank << ":  " << resBuffer[i] << endl;
}
MPI_Gatherv(resBuffer, recvCount[ProcRank], MPI_DOUBLE, r, recvCount, gatherDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(resBuffer, recvCount[ProcRank], MPI_DOUBLE, p, recvCount, gatherDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
*/
/*if (ProcRank == 0)
{
	cout << "r" << endl;
	for (int i = 0; i < size; i++)
		cout << r[i] << "  ";

}*/