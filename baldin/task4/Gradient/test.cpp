#include "test.h"

void GradientTimeTests::SetUp()
{
	workTime = 0;
}

TEST_F(GradientCorrectTests, CheckGradientCorrectSize10)
{
	matr.reset(10, 10);
	matr.fillWithRandomForGradient(rand);

	vec.reset(10);
	vec.fillWithRandom(rand);

	result = gradientMethod(matr, vec, 0.001);
	int check = checkCorrectSolution(matr, vec, result, 0.001);
	
	ASSERT_NE(check, -1); 
}

TEST_F(GradientCorrectTests, CheckGradientCorrectSize100)
{
	matr.reset(100, 100);
	matr.fillWithRandomForGradient(rand);

	vec.reset(100);
	vec.fillWithRandom(rand);

	result = gradientMethod(matr, vec, 0.001);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckGradientCorrectSize500)
{
	matr.reset(500, 500);
	matr.fillWithRandomForGradient(rand);

	vec.reset(500);
	vec.fillWithRandom(rand);

	result = gradientMethod(matr, vec, 0.001);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientTimeTests, TimeTestSize10)
{
	matr.reset(10, 10);
	vec.reset(10);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethod(matr, vec, 0.001);
		endTime = clock();
		workTime += endTime - startTime;
	}
	
	cout << "Average time = " << workTime/10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestSize100)
{
	matr.reset(100, 100);
	vec.reset(100);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethod(matr, vec, 0.001);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestSize500)
{
	matr.reset(500, 500);
	vec.reset(500);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethod(matr, vec, 0.001);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads2Size10)
{
	matr.reset(10, 10);
	matr.fillWithRandomForGradient(rand);

	vec.reset(10);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 2);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads2Size100)
{
	matr.reset(100, 100);
	matr.fillWithRandomForGradient(rand);

	vec.reset(100);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 2);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads2Size500)
{
	matr.reset(500, 500);
	matr.fillWithRandomForGradient(rand);

	vec.reset(500);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 2);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads4Size10)
{
	matr.reset(10, 10);
	matr.fillWithRandomForGradient(rand);

	vec.reset(10);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 4);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads4Size100)
{
	matr.reset(100, 100);
	matr.fillWithRandomForGradient(rand);

	vec.reset(100);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 4);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckOMPGradientCorrectThreads4Size500)
{
	matr.reset(500, 500);
	matr.fillWithRandomForGradient(rand);

	vec.reset(500);
	vec.fillWithRandom(rand);

	result = gradientMethodOMP(matr, vec, 0.001, 4);
	int check = checkCorrectSolution(matr, vec, result, 0.001);

	ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads2Size10)
{
    matr.reset(10, 10);
    matr.fillWithRandomForGradient(rand);

    vec.reset(10);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 2);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads2Size100)
{
    matr.reset(100, 100);
    matr.fillWithRandomForGradient(rand);

    vec.reset(100);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 2);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads2Size500)
{
    matr.reset(500, 500);
    matr.fillWithRandomForGradient(rand);

    vec.reset(500);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 2);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads4Size10)
{
    matr.reset(10, 10);
    matr.fillWithRandomForGradient(rand);

    vec.reset(10);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 4);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads4Size100)
{
    matr.reset(100, 100);
    matr.fillWithRandomForGradient(rand);

    vec.reset(100);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 4);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientCorrectTests, CheckTBBGradientCorrectThreads4Size500)
{
    matr.reset(500, 500);
    matr.fillWithRandomForGradient(rand);

    vec.reset(500);
    vec.fillWithRandom(rand);

    result = gradientMethodTBB(matr, vec, 0.001, 4);
    int check = checkCorrectSolution(matr, vec, result, 0.001);

    ASSERT_NE(check, -1);
}

TEST_F(GradientTimeTests, TimeTestOMPThreads2Size10)
{
	matr.reset(10, 10);
	vec.reset(10);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethodOMP(matr, vec, 0.001, 2);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads2Size100)
{
	matr.reset(100, 100);
	vec.reset(100);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethodOMP(matr, vec, 0.001, 2);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads2Size500)
{
	matr.reset(500, 500);
	vec.reset(500);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethodOMP(matr, vec, 0.001, 2);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads4Size10)
{
	matr.reset(10, 10);
	vec.reset(10);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethodOMP(matr, vec, 0.001, 4);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads4Size100)
{
	matr.reset(100, 100);
	vec.reset(100);
	for (int i = 0; i < 10; ++i)
	{
		matr.fillWithRandomForGradient(rand);
		vec.fillWithRandom(rand);
		startTime = clock();
		result = gradientMethodOMP(matr, vec, 0.001, 4);
		endTime = clock();
		workTime += endTime - startTime;
	}

	cout << "Average time = " << workTime / 10 << " ms." << endl;

	SUCCEED();
}


TEST_F(GradientTimeTests, TimeTestOMPThreads4Size500)
{
    matr.reset(500, 500);
    vec.reset(500);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodOMP(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads2Size1000)
{
    matr.reset(1000, 1000);
    vec.reset(1000);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodOMP(matr, vec, 0.001, 2);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestOMPThreads4Size1000)
{
    matr.reset(1000, 1000);
    vec.reset(1000);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodOMP(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads2Size10)
{
    matr.reset(10, 10);
    vec.reset(10);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 2);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads2Size100)
{
    matr.reset(100, 100);
    vec.reset(100);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 2);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads2Size500)
{
    matr.reset(500, 500);
    vec.reset(500);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 2);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads4Size10)
{
    matr.reset(10, 10);
    vec.reset(10);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads4Size100)
{
    matr.reset(100, 100);
    vec.reset(100);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}


TEST_F(GradientTimeTests, TimeTestTBBThreads4Size500)
{
    matr.reset(500, 500);
    vec.reset(500);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads2Size1000)
{
    matr.reset(1000, 1000);
    vec.reset(1000);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 2);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}

TEST_F(GradientTimeTests, TimeTestTBBThreads4Size1000)
{
    matr.reset(1000, 1000);
    vec.reset(1000);
    for (int i = 0; i < 10; ++i)
    {
        matr.fillWithRandomForGradient(rand);
        vec.fillWithRandom(rand);
        startTime = clock();
        result = gradientMethodTBB(matr, vec, 0.001, 4);
        endTime = clock();
        workTime += endTime - startTime;
    }

    cout << "Average time = " << workTime / 10 << " ms." << endl;

    SUCCEED();
}


