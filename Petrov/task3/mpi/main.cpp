#include <mpi.h>
#include <iostream>
#include <ctime>
using namespace std;

int ProcNum = 0;
int ProcRank = 0;
int GridSize; // размер решетки
int GridCoord[2]; // координыты процессора в решетке
MPI_Comm GridComm; // коммуникатор сетки
MPI_Comm ColComm; // коммуникатор столбцов
MPI_Comm RowComm; // коммуникатор строк

void MatrixInit(double* MatrixA, double* MatrixB, int Size) {
	
		int i, j; // Loop variables
		srand(unsigned(time(0)));
		for (i = 0; i < Size; i++)
			for (j = 0; j < Size; j++) {
				MatrixA[i * Size + j] = rand() / double(1000);
				MatrixB[i * Size + j] = rand() / double(1000);
			}
}

void PrintMatrix(double* Matrix, int Size)
{
	for (int i = 0; i < Size; ++i) {
		for (int j = 0; j < Size; ++j) {
			cout << Matrix[i * Size + j] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}

int CorrectCalc(double*& MatrixA, double*& MatrixB, double*& MatrixC, int Size) {
	if (ProcRank == 0) {
		double temp = 0;
		double* TempMatrix = new double[Size * Size];
		for (int i = 0; i < Size * Size; i++) {
			TempMatrix[i] = 0;
		}
		for (int i = 0; i < Size; i++) {
			for (int j = 0; j < Size; j++) {
				temp = 0;
				for (int k = 0; k < Size; k++)
					temp += MatrixA[i * Size + k] * MatrixB[k * Size + j];
				TempMatrix[i * Size + j] += temp;
			}
		}
		cout << "Comparsion matrix: " << endl;
		PrintMatrix(TempMatrix, Size);
		for (int i = 0; i < Size * Size; i++) {
			if (MatrixC[i] != TempMatrix[i]) {
				delete[]TempMatrix;
				return 0;
			}
		}
		delete[]TempMatrix;
		return 1;
	}
}

// собираем блоки в результат
void ResultCollection(double* MatrixC, double* BlockC, int Size, int BlockSize) {
	double* ResultRow = new double[Size * BlockSize];
	for (int i = 0; i < BlockSize; i++) {
		MPI_Gather(&BlockC[i * BlockSize], BlockSize, MPI_DOUBLE, &ResultRow[i * Size], BlockSize, MPI_DOUBLE, 0, RowComm);
	}
	if (GridCoord[1] == 0) {
		MPI_Gather(ResultRow, BlockSize * Size, MPI_DOUBLE, MatrixC, BlockSize * Size, MPI_DOUBLE, 0, ColComm);
	}
	delete[] ResultRow;
}

// Умножение матричных блоков
void BlockMult(double* BlockA, double* BlockB, double* BlockC, int BlockSize) {
	// вычисление произведения матричных блоков
	for (int i = 0; i < BlockSize; i++) {
		for (int j = 0; j < BlockSize; j++) {
			double temp = 0;
			for (int k = 0; k < BlockSize; k++)
				temp += BlockA[i * BlockSize + k] * BlockB[k * BlockSize + j];
			BlockC[i * BlockSize + j] += temp;
		}
	}
}

// в начале каждой итерации iter алгоритма для каждой строки процессной решетки выбирается процесс, который будет рассылать свой блок матрицы А
void BlockACommunication(int iter, double* BlockA, double* TempBlock, int BlockSize) {
	// определяем ведущий процессор
	int Pivot = (GridCoord[0] + iter) % GridSize;
	// копируем передаваемый блок в буфер
	if (GridCoord[1] == Pivot) {
		for (int i = 0; i < BlockSize * BlockSize; i++)
			BlockA[i] = TempBlock[i];
	}
	// раскидываем блок
	MPI_Bcast(BlockA, BlockSize * BlockSize, MPI_DOUBLE, Pivot, RowComm);
}

// циклический сдвиг
void BlockBCommunication(double* BlockB, int BlockSize) {
	MPI_Status Status;
	int NextProc = GridCoord[0] + 1;
	if (GridCoord[0] == GridSize - 1) NextProc = 0;

	int PrevProc = GridCoord[0] - 1;
	if (GridCoord[0] == 0)
		PrevProc = GridSize - 1;
	MPI_Sendrecv_replace(BlockB, BlockSize * BlockSize, MPI_DOUBLE, NextProc, 0, PrevProc, 0, ColComm, &Status);
}

void CalcParallelResult(double* BlockA, double* TempBlock, double* BlockB, double* BlockC, int BlockSize) {
	for (int iter = 0; iter < GridSize; iter++) {
		/*рассылка блока матрицы А по строке процессорной решетки*/
		BlockACommunication(iter, BlockA, TempBlock, BlockSize);
		// умножение
		BlockMult(BlockA, BlockB, BlockC, BlockSize);
		/*циклический сдвиг блоков матрицы В вдоль столбца процессорной решетки*/
		BlockBCommunication(BlockB, BlockSize);
	}
}

// функция разложения матрицы на сетку
void GridMatrixScatter(double* pMatrix, double* pMatrixBlock, int Size, int BlockSize) {
	double* MatrixRow = new double[BlockSize * Size];
	/*делим матрицу горизонтальными полосами между процессорами нулевого столбца решетки*/
	if (GridCoord[1] == 0) {
		MPI_Scatter(pMatrix, BlockSize * Size, MPI_DOUBLE, MatrixRow, BlockSize * Size, MPI_DOUBLE, 0, ColComm);
	}
	/*распределение каждой строки горизонтальной полосы матрицы вдоль строк процессорной решетки*/
	for (int i = 0; i < BlockSize; i++) {
		MPI_Scatter(&MatrixRow[i * Size], BlockSize, MPI_DOUBLE, &(pMatrixBlock[i * BlockSize]), BlockSize, MPI_DOUBLE, 0, RowComm);
	}
	delete[] MatrixRow;
}

// функция распределения данных между процессорами
void DataDistribution(double* MatrixA, double* MatrixB, double* TempBlock, double* BlockB, int Size, int BlockSize) {
	GridMatrixScatter(MatrixA, TempBlock, Size, BlockSize);
	GridMatrixScatter(MatrixB, BlockB, Size, BlockSize);
}

// функция инициализации памяти и данных
void ProcessInitialization(double*& MatrixA, double*& MatrixB, double*& MatrixC, double*& BlockA, double*& BlockB, double*& BlockC, double*& pTemporaryAblock, int& Size, int& BlockSize) {
		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		BlockSize = Size / GridSize;
		BlockA = new double[BlockSize * BlockSize];
		BlockB = new double[BlockSize * BlockSize];
		BlockC = new double[BlockSize * BlockSize];
		pTemporaryAblock = new double[BlockSize * BlockSize];
		for (int i = 0; i < BlockSize * BlockSize; i++) {
			BlockC[i] = 0;
		}
		if (ProcRank == 0) {
			MatrixA = new double[Size * Size];
			MatrixB = new double[Size * Size];
			MatrixC = new double[Size * Size];
			MatrixInit(MatrixA, MatrixB, Size);
		}
}

/*функция создает коммуникатор в виде двумерной квадратной решетки, определяет координаты каждого процесса в этой решетке*/
void CreateGridCommunicators() {
	int DimSize[2]; // число процессоров в каждом измерении сетки
	int Periodic[2]; // =1, если размерность должна быть динамична
	int Subdims[2]; // =1, если размерность должна быть фиксированна
	DimSize[0] = GridSize;
	DimSize[1] = GridSize;
	Periodic[0] = 0;
	Periodic[1] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
	//определение координат каждого процессора в решетке
	MPI_Cart_coords(GridComm, ProcRank, 2, GridCoord);
	// создание коммуникатора для строк
	Subdims[0] = 0; // фиксация размерности
	Subdims[1] = 1; // наличие данной размерности в подсетке
	MPI_Cart_sub(GridComm, Subdims, &RowComm);
	// создание коммуникатора для столбцов
	Subdims[0] = 1;
	Subdims[1] = 0;
	// создание множества коммуникаторов для каждой строки и каждого столбца решетки в отдельности
	MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

void main(int argc, char* argv[]) {
	double* MatrixA; 
	double* MatrixB; 
	double* MatrixC; 
	int Size = 4; 
	int BlockSize; 
	double* BlockA; 
	double* BlockB; 
	double* BlockC; 
	double* TempBlock;
	double Start, Finish, Duration;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	GridSize = sqrt((double)ProcNum);
	if (ProcNum != GridSize * GridSize) {
		if (ProcRank == 0) {
			printf("Number of processes must be a square \n");
		}
	}
	else {
		if (ProcRank == 0) {
			printf("Parallel matrix multiplication program\n");
		}
		// создание процессорной сетки
		CreateGridCommunicators();
		// выделение памяти и инициализация матриц
		ProcessInitialization(MatrixA, MatrixB, MatrixC, BlockA, BlockB, BlockC, TempBlock, Size, BlockSize);
		if (ProcRank == 0) {
			cout << "Matrix A:" << endl;
			PrintMatrix(MatrixA, Size);
		}

		if (ProcRank == 0) {
			cout << "Matrix B:" << endl;
			PrintMatrix(MatrixB, Size);
		}

		DataDistribution(MatrixA, MatrixB, TempBlock, BlockB, Size, BlockSize);
		CalcParallelResult(BlockA, TempBlock, BlockB, BlockC, BlockSize);
		ResultCollection(MatrixC, BlockC, Size, BlockSize);

		if (ProcRank == 0) {
			cout << "Matrix C:" << endl;
			PrintMatrix(MatrixC, Size);
		}
		if (ProcRank == 0) {
			if (CorrectCalc(MatrixA, MatrixB, MatrixC, Size) == 1)
				cout << "Calculations are correct" << endl;
			else cout << "Calculations are correct" << endl;
		}
		if (ProcRank == 0) {
			delete[]MatrixA;
			delete[]MatrixB;
			delete[]MatrixC;
			delete[]BlockA;
			delete[]BlockB;
			delete[]BlockC;
		}
	}
	MPI_Finalize();
}