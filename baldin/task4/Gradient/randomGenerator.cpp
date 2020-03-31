#include "randomGenerator.h"

randomGenerator::randomGenerator()
{
    random_device device;
    engine.seed(device());
    distribution = uniform_real_distribution<double>(0.0, 10.0);

    trunkate = pow(10.0, 3);
}

randomGenerator::randomGenerator(const int& t)
{
    random_device device;
    engine.seed(device());
    distribution = uniform_real_distribution<double>(0.0, 10.0);
    trunkate = trunkate = pow(10.0, t);
}

randomGenerator::randomGenerator(double min, double max)
{
    if (min > max)
        swap(min, max);
    else if (min == max)
    {
        max = min + 1;
    }

    random_device device;
    engine.seed(device());
    distribution = uniform_real_distribution<double>(min, max);
   
    trunkate = pow(10.0, 3);

}

randomGenerator::randomGenerator(const int& t, double min, double max)
{
    if (min > max)
        swap(min, max);
    else if (min == max)
    {
        max = min + 1;
    }

    random_device device;
    engine.seed(device());
    distribution = uniform_real_distribution<double>(min, max);

    trunkate = pow(10.0, t);
}
void randomGenerator::setDistribution(double min, double max)
{
    if (min > max)
        swap(min, max);
    else if (min == max)
    {
        max = min + 1;
    }

    distribution = uniform_real_distribution<double>(min, max);
}
void randomGenerator::setTrunkate(const int& t)
{
    trunkate = pow(10.0, t);
}
int randomGenerator::getTrunkate()
{
    int result = 0;
    double tempT = trunkate;
    while (tempT != 1)
    {
        tempT /= 10.0;
        result++;
    }
    return result;
}
double randomGenerator::getRandomDouble()
{  
    return distribution(engine);
}

double randomGenerator::getRandomTrunkateDouble()
{
    double a = distribution(engine);
    cout << a << endl;
    cout << trunc(a * trunkate) << endl;
    cout << trunc(a * trunkate) / trunkate << endl;

    return trunc(a * trunkate) / trunkate;
}

