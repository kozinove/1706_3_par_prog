#pragma once

#include <random>
#include <iostream>

using namespace std;

class randomGenerator
{
private:
    mt19937 engine;
    uniform_real_distribution<double> distribution;
    double trunkate;

public:
    randomGenerator();
    randomGenerator(const int& t);
    randomGenerator(double min, double max);
    randomGenerator(const int& t, double min, double max);

    void setDistribution(double min, double max);
    void setTrunkate(const int& t);

    int getTrunkate();
    double getRandomDouble();
    double getRandomTrunkateDouble();

};

