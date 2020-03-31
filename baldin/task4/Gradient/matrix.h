#pragma once

#include <iostream>
#include <iomanip>
#include "randomGenerator.h"
#include <omp.h>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

using namespace std;
using namespace tbb;

class mathVector
{
private:
	int size;

public:
	double* data;

	mathVector();
	mathVector(const int& s);
	mathVector(const mathVector& v);
	~mathVector();

	void reset(const int& s);

	int getSize();

	void printVector();

	void fillWithHands();
	void fillWithRandom(randomGenerator& rand);
	void fillWithRandomRange(randomGenerator& rand, double min, double max);

	mathVector& operator= (const mathVector& v);
};

class matrix
{
private:
	int rowsCount;
	int colsCount;

public:
	double* data;

	matrix();	
	matrix(const int& r, const int& c);
	matrix(const matrix& m);
	~matrix();

	void reset(const int& r, const int& c);

	int getRows();
	int getCols();

	void printMatrix();
	
	void fillWithHands();
	void fillWithRandom(randomGenerator& rand);
	void fillWithRandomRange(randomGenerator& rand, double min, double max);
	void fillWithRandomForGradient(randomGenerator& rand);
	void fillWithRandomRangeForGradient(randomGenerator& rand,double min, double max);

	matrix& operator= (const matrix& m);


};


mathVector operator* (matrix& m, mathVector& v);
mathVector operator- (mathVector& v1, mathVector& v2);
mathVector operator+ (mathVector& v1, mathVector& v2);
double operator* (mathVector& v1, mathVector& v2);
mathVector operator* (mathVector& v, double num);
mathVector operator* (double num, mathVector& v);


mathVector gaussMethod(matrix& matr, mathVector& rightVec, const double& error);
mathVector gradientMethod(matrix& matr, mathVector& rightVec, const double& error);
int checkCorrectSolution(matrix& matr, mathVector& rightVec, mathVector& result, const double& error);

/*
mathVector gradientMethodOMP(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads);
mathVector gradientMethodOMPIter(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads);

void gradientIterationOMP(matrix& matr, mathVector& result, mathVector& prevResult, mathVector& R, mathVector& prevR, mathVector& P, const int& size, double& alfa, double& beta, const int& numThreads);
*/


mathVector MmultVOMP(const matrix& matr, const mathVector& vec, const int& size, const int& numThreads);
double VmultVOMP(const mathVector& vec1, const mathVector& vec2, const int& size, const int& numThreads);
double getAlfa(const mathVector& R, const mathVector& tempVec, const mathVector& P, const int& size, const int& numThreads);
double getBeta(const mathVector& R, const mathVector& prevR, const int& size, const int& numThreads);

int checkCorrectSolutionOMP(matrix& matr, mathVector& rightVec, mathVector& result, const double& error, const int& numThreads);
mathVector gradientMethodOMP(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads);


class VmultVTBB
{
private:
    mathVector *vec1_, *vec2_;
    double res_;
    int size_;

public:
    explicit VmultVTBB(mathVector *vec1, mathVector *vec2, const int& size) : vec1_(vec1), vec2_(vec2), size_(size), res_(0) {}
    VmultVTBB(const VmultVTBB& vv, split) : vec1_(vv.vec1_), vec2_(vv.vec2_), size_(vv.size_), res_(0) {}

    void reset(mathVector *vec1, mathVector *vec2);

    void operator()(const blocked_range<int>& r);

    void join(const VmultVTBB& vv);

    double result();

};

class GetAlfaTBB
{
private:
    mathVector *R_, *vec_, *P_;
    double res1_, res2_;
    int size_;

public:
    explicit GetAlfaTBB() : R_(nullptr), vec_(nullptr), P_(nullptr), size_(0) {}
    explicit GetAlfaTBB(mathVector *R, mathVector *vec, mathVector *P, const int& size) : R_(R), vec_(vec), P_(P), size_(size), res1_(0), res2_(0) {}
    GetAlfaTBB(const GetAlfaTBB& ga, split) : R_(ga.R_), vec_(ga.vec_), P_(ga.P_), size_(ga.size_), res1_(0), res2_(0) {}

    void reset(mathVector *R, mathVector *vec, mathVector *P, const int& size);

    void operator()(const blocked_range<int>& r);

    void join(const GetAlfaTBB& ga);

    double result();

};

class GetBetaTBB
{
private:
    mathVector *R_, *prevR_;
    double res1_, res2_;
    int size_;

public:
    explicit GetBetaTBB() : R_(nullptr), prevR_(nullptr), size_(0) {}
    explicit GetBetaTBB(mathVector *R, mathVector *prevR, const int& size) : R_(R), prevR_(prevR), size_(size), res1_(0), res2_(0) {}
    GetBetaTBB(const GetBetaTBB& gb, split) : R_(gb.R_), prevR_(gb.prevR_), size_(gb.size_), res1_(0), res2_(0) {}

    void reset(mathVector *R, mathVector *prevR, const int& size);

    void operator()(const blocked_range<int>& r);

    void join(const GetBetaTBB& gb);

    double result();

};

class MmultVTBB
{
private:
    matrix *matr_;
    mathVector *vec_;
    mathVector *res_;
    int size_;

public:
    MmultVTBB(matrix *matr, mathVector *vec, mathVector *res, const int& size) : matr_(matr), vec_(vec), res_(res), size_(size) {}

    void operator()(const blocked_range<int>& r) const;
};

int checkCorrectSolutionTBB(matrix& matr, mathVector& rightVec, mathVector& result, const double& error, const int& numThreads);
mathVector gradientMethodTBB(matrix& matr, mathVector& rightVec, const double& error, const int& numThreads);