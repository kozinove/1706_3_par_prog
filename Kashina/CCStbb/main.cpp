#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include <omp.h>
#include "tbb\tbb.h"

using namespace std;
typedef double type;

struct Task//задача для каждого потока
{
	int pointStart;
	int pointEnd;
	int IDtask;
	Task() :pointStart(0), pointEnd(0), IDtask(0) {}
	Task(int a, int b, int c) :pointStart(a), pointEnd(b), IDtask(c) {}
};

class Matrix
{
protected:
	int row = 0;
	int col = 0;
public:
	std::vector<type> vec;
	Matrix()
	{
		vec = vector<double>();
	}

	Matrix(int col, int row)
	{
		this->row = row;
		this->col = col;
		vec = vector<type>(row * col); //по столбцам
	}

	Matrix(int col, int row, type val)
	{
		this->row = row;
		this->col = col;
		vec = vector<type>(row * col, val);
	}

	int getRow() const
	{
		return row;
	}

	int getCol() const
	{
		return col;
	}

	type* operator [] (int j)
	{
		return &vec[j * col];
	}

	void transposition()
	{
		Matrix A(row, col);
		for (int i = 0; i < col; i++)
			for (int j = 0; j < row; j++)
			{
				A[i][j] = (*this)[j][i];
			}
		std::swap(row, col);
		vec = A.vec;
	}

	bool operator == (const Matrix& m)
	{
		bool res = (row == m.row & col == m.col);
		if (res == false)
			return false;

		for (int i = 0; i < vec.size(); i++)
			if (vec[i] != m.vec[i]) return false;
		return true;
	}

	static void readMatrix(Matrix& A, Matrix& B, int N, string name)//для вызова классом(читает матрицу из файла)
	{
		ifstream input(name);
		input >> A.vec[0];//считывает размер, но его затирают далее
		for (int i = 0; i < A.vec.size() + B.vec.size(); i++)
		{
			if (i < A.vec.size())
				input >> A.vec[i];
			else
			{
				for (int j = 0; j < B.vec.size(); j++)
					input >> B.vec[j];
			}
		}
		input.close();
		A.transposition();
		B.transposition();
	}

	static void writeMatrix(Matrix& A, int N)//записывает матрицу в файл
	{
		A.transposition();
		ofstream input("resMatrix.txt");
		for (int i = 0; i < A.vec.size(); i++)
		{
			if (i < A.vec.size())
				input << A.vec[i] << "\t ";
			if ((i + 1) % N == 0)
				input << endl;
		}
		input.close();
	}

	static void randMatrixdata(Matrix& A, Matrix& B, int N)//генерирует матрицы без записи в файл
	{
		for (int i = 0; i < N * N; i++)
		{
			int prob;//вероятность
			prob = rand() % 100;
			if (prob < 10)
				A.vec[i] = (rand() % 10) + (double)(rand() % 10) / 10;//генерирует число и дробную часть
			else
				A.vec[i] = 0;
		}

		for (int i = 0; i < N * N; i++)
		{
			int prob;//вероятность
			prob = rand() % 100;
			if (prob < 10)
				B.vec[i] = (rand() % 10) + (double)(rand() % 10) / 10;
			else
				B.vec[i] = 0;
		}
	}
};

class CCSMatrix
{
protected:
	vector<type> values;
	vector<int> rows;
	vector<int> pointer;
	int N;
public:
	friend class FunctorTBB;
	CCSMatrix() {}

	CCSMatrix(int n) :N(n) {
		pointer.push_back(0);
	}

	CCSMatrix(Matrix& M) :N(M.getRow())
	{
		pointer = vector<int>(N + 1);
		int lastCount = 0;
		for (int j = 0; j < N; j++)
		{
			int curCount = 0;
			for (int i = 0; i < N; i++)
			{
				type el = M[j][i];//берем элемент из поступившей матрицы
				if (el != 0.0)
				{
					curCount++;
					rows.push_back(i);//сохраняем в какой строке
					values.push_back(el);
				}
			}
			if (j != 0)
			{
				pointer[j] = pointer[j - 1] + lastCount;//индекс в values, начала в матрице
			}
			lastCount = curCount;
		}
		pointer[N] = (int)values.size();
	}

	void transposition()
	{
		vector<vector<pair<type, int>>> tmp(N, vector<pair<type, int>>());
		int pCount = 0;//номер столбца
		int numElementsInCol = 0;//сколько ненулевых элементов в столбце
		for (int i = 0; i < values.size(); i += numElementsInCol)
		{
			numElementsInCol = pointer[pCount + 1] - pointer[pCount];//вычисление количества элементов в столбце
			for (int z = i; z < i + numElementsInCol; z++)//цикл для просмотра одного столбца
			{
				int row = rows[z];//в какой строке ненулевой
				int col = pCount;
				type el = values[z];
				tmp[row].push_back(make_pair(el, col));
			}
			pCount++;
		}
		pCount = 0;//для pointer
		vector<int>* cols = &rows;
		int lastCount = 0;
		int vCount = 0;
		for (int i = 0; i < N; i++)
		{
			int numElementInCol = (int)tmp[i].size();//размер строки в матрице
			if (numElementInCol > 0)//если в данной строке есть элементы
			{
				for (int j = 0; j < numElementInCol; j++)//заполнение массивов для CCS
				{
					values[vCount] = tmp[i][j].first;//берем элемент
					(*cols)[vCount] = tmp[i][j].second;
					vCount++;
				}
			}
			if (pCount != 0)
			{
				pointer[pCount] = pointer[pCount - 1] + lastCount;
			}
			lastCount = numElementInCol;
			pCount++;

		}
	}
	void unite(const CCSMatrix& m)//объединение частей в одну матрицу
	{
		for (int i = 0; i < m.values.size(); i++)//добавление values и rows
		{
			values.push_back(m.values[i]);
			rows.push_back(m.rows[i]);
		}
		for (int i = 1; i < m.pointer.size(); i++)//вычисление pointer
		{
			int start = pointer[pointer.size() - 1];//какой pointer был последним
			pointer.push_back(m.pointer[i] - m.pointer[i - 1] + start);
		}
	}

	Matrix CCStoMatrix()
	{
		Matrix res(N, N);

		int Elem = 0;
		int count = 0;
		int resCount = 0;

		for (int i = 0; i < N; i++)
		{
			int numOfElem = pointer[i + 1] - pointer[i];

			if (numOfElem != 0)
				for (int j = 0; j < N; j++)
				{
					if (numOfElem != 0 && rows[resCount] == j)
					{
						res.vec[count] = values[Elem];
						Elem++;
						resCount++;
						numOfElem--;
					}
					else
						res.vec[count] = 0;
					count++;
				}
			else
				for (int j = 0; j < N; j++)
				{
					res.vec[count] = 0;
					count++;
				}
		}

		return res;
	}


	CCSMatrix parallelMult(CCSMatrix& m, const int numThreads)
	{
		class FunctorTBB//Класс-функтор, 
			//принимающий на вход две матрицы и оперирующий рабочими векторами
		{
		private:
			CCSMatrix* m1;
			CCSMatrix* m2;
			vector<Task>* tasks;
			vector<CCSMatrix>* tmp;
			vector<int>* elCountM;
		public:
			FunctorTBB()
			{
				m1 = nullptr;
				m2 = nullptr;
				tasks = nullptr;
				tmp = nullptr;
				elCountM = nullptr;
			}
			FunctorTBB(CCSMatrix* m_1, CCSMatrix* m_2, vector<Task>* tasks_, vector<CCSMatrix>* tmp_, vector<int>* elCountM1)
			{
				m1 = m_1;
				m2 = m_2;
				tasks = tasks_;
				tmp = tmp_;
				elCountM = elCountM1;
			}
			FunctorTBB(const FunctorTBB& f)
			{
				m1 = f.m1;
				m2 = f.m2;
				tasks = f.tasks;
				tmp = f.tmp;
				elCountM = f.elCountM;
			}
			void operator()(const tbb::blocked_range<size_t>& range) const //входом оператора является итерационное пространство
			{
				int N = m1->N;
				vector<int>* cols = &(m1->rows);
				for (size_t th = range.begin(); th < range.end(); th++)//каждый th это поток
					for (int j = (*tasks)[th].pointStart; j < (*tasks)[th].pointEnd; j++)//проходим по ненулевым элементам столбцов
					{
						int numElInResCol = 0;//количество элементов в результирующем столбце
						const int numElCol = m2->pointer[j + 1] - m2->pointer[j];//количество элементов в столбце
						if (numElCol == 0)//если столбец только с нулями
						{
							int size = (*tmp)[th].pointer.size();
							(*tmp)[th].pointer.push_back((*tmp)[th].pointer[size - 1]);
							continue;
						}
						int elCountThis = 0;
						for (int i = 0; i < N; i++)
						{
							const int numElRow = m1->pointer[i + 1] - m1->pointer[i];
							if (numElRow == 0)
							{
								continue;
							}
							int tmpNumElCol = numElCol;//для B
							int tmpNumElRow = numElRow;//для A

							type sum = 0;
							int tmpElCountM = (*elCountM)[th];
							for (int z = 0; z < min(tmpNumElCol, tmpNumElRow);) // проходим по столбцу и умножаем
							{
								int colThis = (*cols)[elCountThis];//cols==rows(из A)
								int rowM = m2->rows[tmpElCountM];//индекс rows с которого нужно начинать task(получаем значение rows в B)
								if (colThis == rowM)
								{
									sum += m1->values[elCountThis] * m2->values[tmpElCountM];
									tmpNumElCol--;
									tmpNumElRow--;
									tmpElCountM++;
									elCountThis++;
								}
								else if (colThis < rowM)
								{
									tmpNumElRow--;//уменьшаем количество элементов
									elCountThis++;//передвигаем индекс
								}
								else
								{
									tmpNumElCol--;
									tmpElCountM++;
								}
							}
							for (int z = 0; z < tmpNumElRow; z++)//двигаем индекс в rows, если у одной матрицы в столбце меньше чисел
								elCountThis++;

							if (sum != 0)//записываем в результирующую матрицу(по столбцам)
							{
								(*tmp)[th].values.push_back(sum);
								(*tmp)[th].rows.push_back(i);
								numElInResCol++;
							}
						}
						const int size = (*tmp)[th].pointer.size();
						(*tmp)[th].pointer.push_back((*tmp)[th].pointer[size - 1] + numElInResCol);
						(*elCountM)[th] += numElCol;
					}
			}
			~FunctorTBB()
			{

			}
		};
		
		int taskNum = numThreads;
		if (taskNum > N)//если потоков больше, чем размерность матрицы
			taskNum = N;
		vector<CCSMatrix> tmp(taskNum, CCSMatrix(N));
		vector<int> elCountM(taskNum);
		vector<int>* cols = &rows;
		vector<Task> task(taskNum);//получаем вектор, в котором taskNum задач
		//создание задач потокам
		{
			int taskSize = m.N / taskNum + (bool)(m.N % taskNum);//если кратное, то + 0, иначе + 1 (распределяем столбцы потокам)
			for (int i = 0; i < taskNum; i++)
				task[i] = Task(i * taskSize, min((i + 1) * taskSize, m.N), i % taskNum);//(начало, конец, индекс)
			int MlastPoint = 0;
			for (int i = 0; i < taskNum; i++)
			{
				elCountM[i] = MlastPoint;
				int starts = task[i].pointStart;
				const int ends = task[i].pointEnd;
				for (starts; starts < ends; starts++)
				{
					MlastPoint += m.pointer[starts + 1] - m.pointer[starts];//считаем количество ненулевых элементов
				}
			}
		}

		FunctorTBB f(this, &m, &task, &tmp, &elCountM);
		tbb::parallel_for(tbb::blocked_range<size_t>((size_t)0, task.size(), (size_t)4), f);//Распараллеливание цикла
		//итерационное пространство blocked_range
		for (int i = 1; i < tmp.size(); i++)
		{
			tmp[0].unite(tmp[i]);
		}
		if (tmp[0].pointer.size() < N + 1)
			tmp[0].pointer.push_back(tmp[0].values.size());
		return tmp[0];
	}
	   	
};


int main(int argc, char* argv[])
{
	int numThreads;
	int mSize;
	cout << "Number of Threads = ";
	cin >> numThreads;

	tbb::task_scheduler_init init(numThreads); //Инициализация библиотеки

	cout << "1.Random matrix" << endl;
	cout << "2.Read from file" << endl;
	int choice;
	cin >> choice;

	if (choice == 1)
	{
		cout << "Size of Matrix = ";
		cin >> mSize;
		Matrix a(mSize, mSize), b(mSize, mSize), c(mSize, mSize);//создаем матрицы
		Matrix::randMatrixdata(a, b, mSize);//сгенерировали матрицы
		CCSMatrix a1(a), b1(b), c1(c);//преобразовали обычную в CCS матрицы

		a1.transposition();

		double time = omp_get_wtime();//подсчет времени
		c1 = a1.parallelMult(b1, numThreads);
		time = omp_get_wtime() - time;

		cout << time << endl;//завершение подсчета времени и вывод

		c = c1.CCStoMatrix();
		Matrix::writeMatrix(c, mSize);//запись результата в файл

	}
	else if (choice == 2)
	{
		ifstream input("matrix.txt");//читаем из файла размерность матрицы
		input >> mSize;
		input.close();

		Matrix a(mSize, mSize), b(mSize, mSize), c(mSize, mSize);//создаем матрицы
		Matrix::readMatrix(a, b, mSize, "matrix.txt");
		CCSMatrix a1(a), b1(b), c1(c);

		a1.transposition();

		double time = omp_get_wtime();//подсчет времени
		c1 = a1.parallelMult(b1, numThreads);
		time = omp_get_wtime() - time;

		cout << time << endl;//завершение подсчета времени и вывод

		c = c1.CCStoMatrix();
		Matrix::writeMatrix(c, mSize);//запись результата в файл
	}
	else
		cout << "Error" << endl;

	return 0;
}