#include "matrix.h"

matrix::matrix()
{
    rowsCount = 1;
    colsCount = 1;

    data = new double[1];
    data[0] = 0;
}

matrix::matrix(const int& r, const int& c)
{
    rowsCount = r;
    colsCount = c;

    data = new double[rowsCount * colsCount];
    for (int i = 0; i < rowsCount * colsCount; ++i)
    {
        data[i] = 0;
    }
}

matrix::matrix(const matrix& m)
{
    rowsCount = m.rowsCount;
    colsCount = m.colsCount;

    data = new double[rowsCount * colsCount];

    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            data[i * colsCount + j] = m.data[i * colsCount + j];
}

matrix::~matrix()
{
    delete[] data;
}

void matrix::reset(const int& r, const int& c)
{
    rowsCount = r;
    colsCount = c;

    delete[] data;

    data = new double[rowsCount * colsCount];

    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            data[i * colsCount + j] = 0;
}

int matrix::getRows()
{
    return rowsCount;
}

int matrix::getCols()
{
    return colsCount;
}

void matrix::printMatrix()
{
    cout << "Matrix: " << endl;
    for (int i = 0; i < rowsCount; ++i)
    {
        for (int j = 0; j < colsCount; ++j)
        {
            cout << data[i * colsCount + j] << "   ";
        }
        cout << endl;
    }
}

void matrix::fillWithHands()
{
    cout << "Write matrix: " << endl;
    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            cin >> data[i * colsCount + j];
}

void matrix::fillWithRandom(randomGenerator& rand)
{
    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            data[i * colsCount + j] = rand.getRandomDouble();
}

void matrix::fillWithRandomRange(randomGenerator& rand, double min, double max)
{
    rand.setDistribution(min, max);

    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            data[i * colsCount + j] = rand.getRandomDouble();
}
void matrix::fillWithRandomForGradient(randomGenerator& rand)
{
    if (rowsCount != colsCount)
    {
        cout << "Rows count != cols count. Generating default random matrix." << endl;

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                data[i * colsCount + j] = rand.getRandomDouble();
    }
    else
    {
        double* tempMatrix = new double[rowsCount * colsCount];

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                tempMatrix[i * colsCount + j] = rand.getRandomDouble();

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                data[i * colsCount + j] = sqrt(tempMatrix[i * colsCount + j] * tempMatrix[j * rowsCount + i]);

        delete[] tempMatrix;
    }
}
void matrix::fillWithRandomRangeForGradient(randomGenerator& rand, double min, double max)
{
    rand.setDistribution(min, max);

    if (rowsCount != colsCount)
    {
        cout << "Rows count != cols count. Generating default random matrix." << endl;

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                data[i * colsCount + j] = rand.getRandomDouble();
    }
    else
    {
        double* tempMatrix = new double[rowsCount * colsCount];

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                tempMatrix[i * colsCount + j] = rand.getRandomDouble();

        for (int i = 0; i < rowsCount; ++i)
            for (int j = 0; j < colsCount; ++j)
                data[i * colsCount + j] = sqrt(tempMatrix[i * colsCount + j] * tempMatrix[j * colsCount + i]);

        delete[] tempMatrix;
    }
}

matrix& matrix::operator=(const matrix& m)
{
    rowsCount = m.rowsCount;
    colsCount = m.colsCount;
    delete[] data;
    data = new double[rowsCount * colsCount];
    for (int i = 0; i < rowsCount; ++i)
        for (int j = 0; j < colsCount; ++j)
            data[i * colsCount + j] = m.data[i * colsCount + j];
    return *this;
}

mathVector::mathVector()
{
    size = 1;

    data = new double[1];
    data[0] = 0;
}

mathVector::mathVector(const int& s)
{
    size = s;

    data = new double[size];

    for (int i = 0; i < size; ++i)
        data[i] = 0;
}

mathVector::mathVector(const mathVector& v)
{
    size = v.size;
    data = new double[size];
    for (int i = 0; i < size; ++i)
        data[i] = v.data[i];
}

mathVector::~mathVector()
{
    delete[] data;
}

void mathVector::reset(const int& s)
{
    size = s;

    delete[] data;

    data = new double[size];

    for (int i = 0; i < size; ++i)
        data[i] = 0;
}

int mathVector::getSize()
{
    return size;
}

void mathVector::printVector()
{
    cout << "Vector: " << endl;
    for (int i = 0; i < size; ++i)
    {
        cout << data[i] << "   ";
    }
    cout << endl;
}

void mathVector::fillWithHands()
{
    cout << "Write vector: " << endl;
    for (int i = 0; i < size; ++i)
        cin >> data[i];
}

void mathVector::fillWithRandom(randomGenerator& rand)
{
    for (int i = 0; i < size; ++i)
        data[i] = rand.getRandomDouble();
}

void mathVector::fillWithRandomRange(randomGenerator& rand, double min, double max)
{
    rand.setDistribution(min, max);
    for (int i = 0; i < size; ++i)
        data[i] = rand.getRandomDouble();
}

mathVector& mathVector::operator=(const mathVector& v)
{
    size = v.size;
    delete[] data;
    data = new double[size];
    for (int i = 0; i < size; ++i)
        data[i] = v.data[i];
    return *this;
}

mathVector operator* (matrix& m, mathVector& v)
{
    int rows = m.getRows();
    int cols = m.getCols();

    mathVector result(rows);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result.data[i] += m.data[cols * i + j] * v.data[j];
    return result;
}

mathVector operator- (mathVector& v1, mathVector& v2)
{
    int size = v1.getSize();
    if (size != v2.getSize())
    {
        cout << "Vector 1 size != vector 2 size." << endl;
        return mathVector();
    }
    else
    {
        mathVector result(size);

        for (int i = 0; i < size; ++i)
            result.data[i] = v1.data[i] - v2.data[i];

        return result;
    }
}

mathVector operator+ (mathVector& v1, mathVector& v2)
{
    int size = v1.getSize();
    if (size != v2.getSize())
    {
        cout << "Vector 1 size != vector 2 size." << endl;
        return mathVector();
    }
    else
    {
        mathVector result(size);

        for (int i = 0; i < size; ++i)
            result.data[i] = v1.data[i] + v2.data[i];

        return result;
    }
}

double operator* (mathVector& v1, mathVector& v2)
{
    int size = v1.getSize();
    if (size != v2.getSize())
    {
        cout << "Vector 1 size != vector 2 size." << endl;
        return 0;
    }
    else
    {
        double result = 0;

        for (int i = 0; i < size; ++i)
            result += v1.data[i] * v2.data[i];

        return result;
    }
}

mathVector operator* (mathVector& v, double num)
{
    int size = v.getSize();
    mathVector result(size);

    for (int i = 0; i < size; ++i)
        result.data[i] = v.data[i] * num;

    return result;
}

mathVector operator* (double num, mathVector& v)
{
    return v * num;
}

mathVector gaussMethod(matrix& matr, mathVector& rightVec, const double& error)
{

    if (matr.getCols() != matr.getRows() && matr.getCols() != rightVec.getSize())
    {
        cout << "Wrong matrix or vector size, return zero-vector." << endl;
        return mathVector(matr.getRows());
    }

    int size = matr.getCols();
    matrix tempMatr = matr;
    mathVector tempVec = rightVec;
    mathVector result(size);
    double coeff = 0;


    for (int i = 0; i < size; ++i)
    {
        for (int j = i + 1; j < size; ++j)
        {
            coeff = tempMatr.data[j * size + i] / tempMatr.data[i * size + i];
            for (int k = i; k < size; ++k)
            {
                tempMatr.data[j * size + k] -= coeff * tempMatr.data[i * size + k];
            }
            tempVec.data[j] -= coeff * tempVec.data[i];
        }
    }


    result.data[size - 1] = tempVec.data[size - 1] / tempMatr.data[size * size - 1];
    for (int i = size - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < size; ++j)
        {
            sum += tempMatr.data[i * size + j] * result.data[j];
        }
        result.data[i] = (tempVec.data[i] - sum) / tempMatr.data[i * size + i];
    }

    int check;
    check = checkCorrectSolution(matr, rightVec, result, error);
    if (check == 1)
        cout << "Gauss solution is full correct." << endl;
    else if (check == 0)
        cout << "Gauss sulution correct with allowable error." << endl;
    else
        cout << "Wrong gauss solution." << endl;


    return result;
}

mathVector gradientMethod(matrix& matr, mathVector& rightVec, const double& error)
{
    if (matr.getCols() != matr.getRows() && matr.getCols() != rightVec.getSize())
    {
        cout << "Wrong matrix or vector size for Gradient method, return zero-vector." << endl;
        return mathVector(matr.getRows());
    }

    int size = matr.getCols();
    int iteration = 0;

    mathVector Result(size);
    mathVector prevResult(size);
    mathVector R(size);
    mathVector prevR(size);
    mathVector P(size);
    double alfa = 0;
    double beta = 0;



    R = matr * prevResult;
    R = rightVec - R;
    P = R;

    mathVector tempVec1(size);
    mathVector tempVec2(size);
    double tempD1 = 0;
    double tempD2 = 0;

    int addIterations = 0;
    if (size <= 10)
        addIterations = 1;
    else if (size >= 10 && size < 50)
        addIterations = 10;
    else if (size >= 50 && size < 100)
        addIterations = 30;
    else if (size >= 100 && size < 250)
        addIterations = 100;
    else if (size >= 250 && size < 500)
        addIterations = 250;
    else if (size >= 500 && size < 1000)
        addIterations = 500;
    else
        addIterations = size * 2;


    while (iteration < size + addIterations)
    {
        tempVec1 = matr * P;
        alfa = (R * R) / (tempVec1 * P);

        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;

        beta = (R * R) / (prevR * prevR);

        tempVec1 = P * beta;
        P = R + tempVec1;

        //P.printVector();

        iteration++;

        prevResult = Result;


    }

    int check = 0;
    while (true)
    {
        check = checkCorrectSolution(matr, rightVec, Result, error);
        if (check == 0 || check == 1)
            break;

        tempVec1 = matr * P;
        alfa = (R * R) / (tempVec1 * P);

        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;

        beta = (R * R) / (prevR * prevR);

        tempVec1 = P * beta;
        P = R + tempVec1;

        iteration++;

        prevResult = Result;

        if (iteration > size * 10)
        {
            break;
        }
    }

    //cout << "Iterations count: " << iteration << endl;

    /*if (check == 1)
        cout << "Gradient solution is full correct." << endl;
    else if (check == 0)
        cout << "Gradient sulution correct with allowable error." << endl;
    else
        cout << "Wrong gradient solution." << endl;*/


    return Result;
}

int checkCorrectSolution(matrix& matr, mathVector& rightVec, mathVector& result, const double& error)
{
    double currError = 0;

    mathVector residual = rightVec;

    mathVector tempVec = matr * result;

    residual = residual - tempVec;

    for (int i = 0; i < residual.getSize(); ++i)
    {
        if (residual.data[i] != 0)
        {
            currError += abs(residual.data[i]);
        }
    }
    if (currError == 0)
        return 1;
    else if (currError < error)
        return 0;
    else
        return -1;
}


int checkCorrectSolutionOMP(matrix& matr, mathVector& rightVec, mathVector& result, const double& error, const int& numThreads)
{
    double currError = 0;

    mathVector residual = rightVec;

    mathVector tempVec = MmultVOMP(matr, result, rightVec.getSize(), numThreads);

    residual = residual - tempVec;

    for (int i = 0; i < residual.getSize(); ++i)
    {
        if (residual.data[i] != 0)
        {
            currError += abs(residual.data[i]);
        }
    }
    if (currError == 0)
        return 1;
    else if (currError < error)
        return 0;
    else
        return -1;
}


mathVector MmultVOMP(const matrix& matr, const mathVector& vec, const int& size, const int& numThreads)
{
    mathVector result(size);
    int i = 0;
    int j = 0;

#pragma omp parallel for private(j) num_threads(numThreads)
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            result.data[i] += matr.data[size * i + j] * vec.data[j];
        }
    }

    return result;
}

double VmultVOMP(const mathVector& vec1, const mathVector& vec2, const int& size, const int& numThreads)
{
    double result = 0;
    int i = 0;
#pragma omp parallel for private(i) reduction(+:result) num_threads(numThreads)
    for (i = 0; i < size; ++i)
    {
        result += vec1.data[i] * vec2.data[i];
    }

    return result;
}

double getAlfa(const mathVector& R, const mathVector& tempVec, const mathVector& P, const int& size, const int& numThreads)
{
    double result = 0;
    double temp = 0;
    int i = 0;

#pragma omp parallel for private(i) reduction(+:result, temp) num_threads(numThreads)
    for (i = 0; i < size; ++i)
    {
        result += R.data[i] * R.data[i];
        temp += tempVec.data[i] * P.data[i];
    }

    return result / temp;
}

double getBeta(const mathVector& R, const mathVector& prevR, const int& size, const int& numThreads)
{
    double result = 0;
    double temp = 0;
    int i = 0;

#pragma omp parallel for private(i) reduction(+:result, temp) num_threads(numThreads)
    for (i = 0; i < size; ++i)
    {
        result += R.data[i] * R.data[i];
        temp += prevR.data[i] * prevR.data[i];
    }

    return result / temp;
}


mathVector gradientMethodOMP(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads)
{
    if (matr.getCols() != matr.getRows() && matr.getCols() != rightVec.getSize())
    {
        cout << "Wrong matrix or vector size for Gradient method, return zero-vector." << endl;
        return mathVector(matr.getRows());
    }

    int size = matr.getCols();


    mathVector Result(size);
    mathVector prevResult(size);
    mathVector R(size);
    mathVector prevR(size);
    mathVector P(size);
    double alfa = 0;
    double beta = 0;


    int i = 0;
    int j = 0;
    mathVector tempVec1(size);
    mathVector tempVec2(size);
    mathVector residual = rightVec;
    double tempD1 = 0;
    double tempD2 = 0;
    double currError = 0;
    int iteration = 0;


#pragma omp parallel num_threads(numThreads) private(i, j)
    {
#pragma omp for
        for (i = 0; i < size; ++i)
        {
            for (j = 0; j < size; ++j)
            {
                tempVec1.data[i] += matr.data[size * i + j] * prevResult.data[j];
            }
            R.data[i] = rightVec.data[i] - tempVec1.data[i];
            P.data[i] = R.data[i];
        }
    }
   
   

    int addIterations = 0;
    if (size <= 10)
        addIterations = 1;
    else if (size >= 10 && size < 50)
        addIterations = 10;
    else if (size >= 50 && size < 100)
        addIterations = 30;
    else if (size >= 100 && size < 250)
        addIterations = 125;
    else if (size >= 250 && size < 500)
        addIterations = 250;
    else if (size >= 500 && size < 1000)
        addIterations = 500;
    else
        addIterations = size * 2;

    while (iteration < size + addIterations)
    {
        tempVec1 = MmultVOMP(matr, P, size, numThreads);
        //tempVec1.printVector();
        alfa = getAlfa(R, tempVec1, P, size, numThreads);
        //cout << alfa << endl;

        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        /*cout << iteration << endl;
        Result.printVector();*/
        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;

        /*cout << iteration << endl;
        R.printVector();*/

        beta = getBeta(R, prevR, size, numThreads);

        /*cout << iteration << endl;
        cout << beta << endl;*/

        tempVec1 = P * beta;
        P = R + tempVec1;

        /*cout << iteration << endl;
        P.printVector();*/

        prevResult = Result;

        iteration++;
    }

    int check = 0;
    while (true)
    {
        check = checkCorrectSolutionOMP(matr, rightVec, Result, error, numThreads);
        if (check == 0 || check == 1)
            break;

        tempVec1 = MmultVOMP(matr, P, size, numThreads);
        alfa = getAlfa(R, tempVec1, P, size, numThreads);

        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;

        beta = getBeta(R, prevR, size, numThreads);

        tempVec1 = P * beta;
        P = R + tempVec1;

        iteration++;

        prevResult = Result;

        if (iteration > size * 10)
        {
            break;
        }
    }
    return Result;
}

void VmultVTBB::reset(mathVector *vec1, mathVector *vec2)
{
    vec1_ = vec1;
    vec2_ = vec2;
    res_ = 0;
}

void VmultVTBB::operator()(const blocked_range<int>& r)
{
    int begin = r.begin(), end = r.end();
    for (int i = begin; i != end; i++)
        res_ += vec1_->data[i] * vec2_->data[i];
}

void VmultVTBB::join(const VmultVTBB & vv)
{
    res_ += vv.res_;
}

double VmultVTBB::result()
{
    return res_;
}

void GetAlfaTBB::reset(mathVector *R, mathVector *vec, mathVector *P, const int& size)
{
    R_ = R;
    vec_ = vec;
    P_ = P;
    size_ = size;
    res1_ = 0;
    res2_ = 0;
}

void GetAlfaTBB::operator()(const blocked_range<int>& r)
{
    int begin = r.begin(), end = r.end();
    for (int i = begin; i != end; i++)
    {
        res1_ += R_->data[i] * R_->data[i];
        res2_ += vec_->data[i] * P_->data[i];
    }
}

void GetAlfaTBB::join(const GetAlfaTBB & ga)
{
    res1_ += ga.res1_;
    res2_ += ga.res2_;
}

double GetAlfaTBB::result()
{
    return res1_ / res2_;
}

void GetBetaTBB::reset(mathVector *R, mathVector *prevR, const int& size)
{
    R_ = R;
    prevR_ = prevR;
    size_ = size;
    res1_ = 0;
    res2_ = 0;
}

void GetBetaTBB::operator()(const blocked_range<int>& r)
{
    int begin = r.begin(), end = r.end();
    for (int i = begin; i != end; i++)
    {
        res1_ += R_->data[i] * R_->data[i];
        res2_ += prevR_->data[i] * prevR_->data[i];
    }
}

void GetBetaTBB::join(const GetBetaTBB & gb)
{
    res1_ += gb.res1_;
    res2_ += gb.res2_;
}

double GetBetaTBB::result()
{
    return res1_ / res2_;
}

void MmultVTBB::operator()(const blocked_range<int>& r) const
{
    int begin = r.begin(), end = r.end();
    for (int i = begin; i != end; i++)
    {
        res_->data[i] = 0;
        for (int j = 0; j < size_; j++)
        {
            res_->data[i] += matr_->data[size_ * i + j] * vec_->data[j];
        }
    }
}

int checkCorrectSolutionTBB(matrix& matr, mathVector& rightVec, mathVector& result, const double& error, const int& numThreads)
{
    double currError = 0;

    mathVector residual = rightVec;

    mathVector tempVec(residual.getSize());

    parallel_for(blocked_range<int>(0, residual.getSize()), MmultVTBB(&matr, &result, &tempVec, residual.getSize()));

    residual = residual - tempVec;

    for (int i = 0; i < residual.getSize(); ++i)
    {
        if (residual.data[i] != 0)
        {
            currError += abs(residual.data[i]);
        }
    }
    if (currError == 0)
        return 1;
    else if (currError < error)
        return 0;
    else
        return -1;
}

mathVector gradientMethodTBB(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads)
{
    if (matr.getCols() != matr.getRows() && matr.getCols() != rightVec.getSize())
    {
        cout << "Wrong matrix or vector size for Gradient method, return zero-vector." << endl;
        return mathVector(matr.getRows());
    }

    task_scheduler_init init(numThreads);

    int size = matr.getCols();

    mathVector Result(size);
    mathVector prevResult(size);
    mathVector R(size);
    mathVector prevR(size);
    mathVector P(size);
    double alfa = 0;
    double beta = 0;


    int i = 0;
    int j = 0;
    mathVector tempVec1(size);
    mathVector tempVec2(size);
    mathVector residual = rightVec;
    double tempD1 = 0;
    double tempD2 = 0;
    double currError = 0;
    int iteration = 0;

    
    parallel_for(blocked_range<int>(0, size), MmultVTBB(&matr, &prevResult, &tempVec1, size));
    R = rightVec - tempVec1;
    P = R;
   

    int addIterations = 0;
    if (size <= 10)
        addIterations = 1;
    else if (size >= 10 && size < 50)
        addIterations = 10;
    else if (size >= 50 && size < 100)
        addIterations = 30;
    else if (size >= 100 && size < 250)
        addIterations = 100;
    else if (size >= 250 && size < 500)
        addIterations = 250;
    else if (size >= 500 && size < 1000)
        addIterations = 500;
    else
        addIterations = size * 2;



    GetAlfaTBB tempGetAlfaTBB;
    GetBetaTBB tempGetBetaTBB;

    while (iteration < size + addIterations)
    {
        parallel_for(blocked_range<int>(0, size), MmultVTBB(&matr, &P, &tempVec1, size));

        tempGetAlfaTBB.reset(&R, &tempVec1, &P, size);
        parallel_reduce(blocked_range<int>(0, size), tempGetAlfaTBB);
        alfa = tempGetAlfaTBB.result();


        //cout << alfa << endl;


        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        /*cout << iteration << endl;
        Result.printVector();*/

        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;
 

        tempGetBetaTBB.reset(&R, &prevR, size);
        parallel_reduce(blocked_range<int>(0, size), tempGetBetaTBB);
        beta = tempGetBetaTBB.result();


        tempVec1 = P * beta;
        P = R + tempVec1;

        /*cout << iteration << endl;
        P.printVector();*/

        prevResult = Result;

        iteration++;

    }

    int check = 0;
    while (true)
    {
        check = checkCorrectSolutionTBB(matr, rightVec, Result, error, numThreads);
        if (check == 0 || check == 1)
            break;

        parallel_for(blocked_range<int>(0, size), MmultVTBB(&matr, &P, &tempVec1, size));

        tempGetAlfaTBB.reset(&R, &tempVec1, &P, size);
        parallel_reduce(blocked_range<int>(0, size), tempGetAlfaTBB);
        alfa = tempGetAlfaTBB.result();

        tempVec2 = P * alfa;
        Result = prevResult + tempVec2;

        prevR = R;

        tempVec1 = tempVec1 * alfa;
        R = R - tempVec1;

        tempGetBetaTBB.reset(&R, &prevR, size);
        parallel_reduce(blocked_range<int>(0, size), tempGetBetaTBB);
        beta = tempGetBetaTBB.result();

        tempVec1 = P * beta;
        P = R + tempVec1;


        prevResult = Result;

        iteration++;

        if (iteration > size * 10)
        {
            break;
        }
    }

    return Result;
}
