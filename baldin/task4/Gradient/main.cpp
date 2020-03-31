#include "matrix.h"
#include "randomGenerator.h"
#include "ui.h"
#include <iostream>
#include <gtest/gtest.h>


using namespace std;

void myTests()
{
    cout << "TBB Tests" << endl;

    task_scheduler_init init(4);
    randomGenerator rand;

    matrix a(1000, 1000);
    a.fillWithRandomForGradient(rand);
   // a.printMatrix();
    mathVector b(1000);
    b.fillWithRandom(rand);
    mathVector c(1000);
    c.fillWithRandom(rand);
    mathVector d(1000);
    d.fillWithRandom(rand);
    //b.printVector();

    cout << "--------------" << endl;
    
    mathVector Q = gradientMethodTBB(a, b, 0.001, 2);
    mathVector W = gradientMethodOMP(a, b, 0.001, 2);

    cout << checkCorrectSolutionTBB(a, b, Q, 0.005, 2) << endl;
    cout << checkCorrectSolutionOMP(a, b, W, 0.005, 2) << endl;



    cout << "Finish" << endl; 
}

int main(int argc, char* argv[])
{
    //cout << fixed << setprecision(30);
    string arg = "";
    if (argc > 1)
        arg = argv[1];

    if (argc == 2 && arg == "-ui")
    {
        randomGenerator rand;
        gradientUI(rand);
        return 1;
    }
    else if (argc == 2 && arg == "-test")
    {
        myTests();
    }
    else if (argc > 1 && arg == "-gtest")
    {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
    else if (argc > 1 && arg == "-data")
    {
        cout << argc << endl;
        randomGenerator rand;
        int size = 10;
        double min = 0.0;
        double max = 10.0;
        double error = 0.0001;
        if (argc > 2)
            size = atoi(argv[2]);
        if (argc > 3)
            min = atof(argv[3]);
        if (argc > 4)
            max = atof(argv[4]);
        if (argc > 5)
            error = atof(argv[5]);


        cout << "Using data: " << endl;
        cout << "Size: " << size << endl;
        cout << "Minimum: " << min << endl;
        cout << "Maximim: " << max << endl;
        cout << "Allowable error: " << error << endl;



        matrix matr(size, size);
        matr.fillWithRandomRangeForGradient(rand, min, max);

        mathVector vec(size);
        vec.fillWithRandomRange(rand, min, max);

        mathVector result = gradientMethod(matr, vec, error);
        return 2;
    }
    else if (argc > 1 && arg == "-help")
    {
        cout << "Available arguments: " << endl;
        cout << "'-help' - Argument list." << endl;
        cout << "'-ui' - Use user interface." << endl;
        cout << "'-test' - Run user tests." << endl;
        cout << "'-gtest' - Run Google tests." << endl;
        cout << "'-data' 'size' 'min' 'max' 'error' - Run Gradient method with written data." << endl;
    }
    else
    {
        if (argc == 1)
        {
            cout << "Use the arguments." << endl;
            cout << "Available arguments: " << endl;
            cout << "'-help' - Argument list." << endl;
            cout << "'-ui' - Use user interface." << endl;
            cout << "'-test' - Run user tests." << endl;
            cout << "'-gtest' - Run Google tests." << endl;
            cout << "'-data' 'size' 'min' 'max' 'error' - Run Gradient method with written data." << endl;
            return -1;
        }
        else
        {
            cout << "Unknown arguments." << endl;
            cout << "Available arguments: " << endl;
            cout << "'-help' - Argument list." << endl;
            cout << "'-ui' - Use user interface." << endl;
            cout << "'-test' - Run user tests." << endl;
            cout << "'-gtest' - Run Google tests." << endl;
            cout << "'-data' 'size' 'min' 'max' 'error' - Run Gradient method with written data." << endl;
            return -2;
        }
    }





    /*double a = 0;
    a = rand.getRandomTrunkateDouble();
    cout << a << endl;*/


    return 0;
}