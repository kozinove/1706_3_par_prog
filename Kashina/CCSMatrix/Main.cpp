#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

typedef double type;

using namespace std;

class Matrix
{
protected:
	int row = 0;
	int col = 0;
public:
	vector<type> vec;
	Matrix()
	{
		vec = vector<double>();
	}
	Matrix(int col, int row)
	{
		this->row = row;
		this->col = col;
		vec = vector<type>(row * col);//по столбцам
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
	type* getP()
	{
		return &vec[0];
	}
	type* operator [] (int j)
	{
		return &vec[j * col];
	}
	Matrix operator * (Matrix& m)
	{
		if (getCol() != m.getRow()) return Matrix();
		Matrix newMatr(getRow(), m.getCol());
		for (int j = 0; j < m.getRow(); j++)
			for (int i = 0; i < getCol(); i++)
				for (int z = 0; z < getCol(); z++)
					newMatr[j][i] += (*this)[z][i] * m[j][z];
		return newMatr;
	}
	void transposition()
	{
		Matrix A(row, col);
		for (int i = 0; i < col; i++)
			for (int j = 0; j < row; j++)
			{
				A[i][j] = (*this)[j][i];
			}
		swap(row, col);
		vec = A.vec;
	}

	bool operator == (const Matrix& m)
	{
		bool res = ((row == m.row) & (col == m.col));
		if (res == false) return false;
		for (int i = 0; i < vec.size(); i++)
			if (vec[i] != m.vec[i]) return false;
		return true;
	}

	static void readMatrix(Matrix& A, Matrix& B, int N)
	{
		ifstream input("source.txt");
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


	static void writeMatrix(Matrix& A, int N)
	{
		A.transposition();
		ofstream input("result.txt");
		for (int i = 0; i < A.vec.size(); i++)
		{
			if (i < A.vec.size())
				input << A.vec[i] << "\t ";
			if ((i + 1) % N == 0)
				input << endl;
		}
		input.close();
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

	CCSMatrix(int n) :N(n) {}
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
	CCSMatrix operator * (const CCSMatrix& m)
	{
		CCSMatrix res(N);//результирующая матрица
		res.pointer.push_back(0);
		transposition();
		vector<int>* cols = &rows;
		int elCountM = 0;
		for (int j = 0; j < m.N; j++)
		{
			int numElInResCol = 0;
			const int numElementInCol = m.pointer[j + 1] - m.pointer[j];//вычисляется количество элементов в столбце
			if (numElementInCol == 0)
			{
				int size = res.pointer.size();//размер pointer в результирующем векторе
				res.pointer.push_back(res.pointer[size - 1]);//дублируем предыдущий(если в данном столбце ноль)
				continue;
			}
			int elCountThis = 0;
			for (int i = 0; i < N; i++)
			{
				const int numElementInRow = pointer[i + 1] - pointer[i];//сколько элементов в строке
				if (numElementInRow == 0)//если ноль элементов
				{
					continue;
				}
				int tmpNumElCol = numElementInCol;//количество элементов в столбце матрицы m
				int tmpNumElRow = numElementInRow;//количество элементов в строке матрицы this

				type sum = 0;
				int tmpElCountM = elCountM;
				for (int z = 0; z < min(tmpNumElCol, tmpNumElRow);)//выбираем min, т.к. соответствующие элементы умножаются, в других будет ноль
				{
					int colThis = (*cols)[elCountThis];//получает позицию ненулевого
					int rowM = m.rows[tmpElCountM];//получает позицию ненулевого в M
					if (colThis == rowM)//сравнивается на соответсвие
					{
						sum += values[elCountThis] * m.values[tmpElCountM];
						tmpNumElCol--;
						tmpNumElRow--;
						tmpElCountM++;
						elCountThis++;
					}
					else if (colThis < rowM)//this
					{
						tmpNumElRow--;//количество оставшихся ненулевых
						elCountThis++;//позиция
					}
					else//для m
					{
						tmpNumElCol--;
						tmpElCountM++;
					}
				}
				for (int z = 0; z < tmpNumElRow; z++)
					elCountThis++;

				if (sum != 0)
				{
					res.values.push_back(sum);
					res.rows.push_back(i);
					numElInResCol++;//сколько элементов в столбце
				}
			}
			const int size = res.pointer.size();
			res.pointer.push_back(res.pointer[size - 1] + numElInResCol);
			elCountM += numElementInCol;
		}
		transposition();
		return res;
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

};



int main()
{
	Matrix a(4, 4), b(4, 4), c(4, 4);

	Matrix::readMatrix(a, b, 4);

	CCSMatrix a1(a), b1(b), c1(c);
	c1 = a1 * b1;
	c = c1.CCStoMatrix();

	Matrix::writeMatrix(c, 4);

	return 0;
}