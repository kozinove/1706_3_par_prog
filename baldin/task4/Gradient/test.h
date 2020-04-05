#pragma once

#include "matrix.h"
#include "randomGenerator.h"
#include <iostream>
#include <ctime> 
#include <gtest/gtest.h>
#include <omp.h>


using namespace std;

class GradientCorrectTests : public ::testing::Test 
{
public:
	randomGenerator rand;
	matrix matr;
	mathVector vec;
	mathVector result;
};

class GradientTimeTests : public ::testing::Test
{
public:
	randomGenerator rand;
	matrix matr;
	mathVector vec;
	mathVector result;

	unsigned int startTime;
	unsigned int endTime;
	double workTime;

	void SetUp();
};
