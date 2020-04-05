#include "ui.h"

void gradientUI(randomGenerator& rand)
{
	int progRun = 1;
	while (progRun)
	{
		cout << "Enter matrix size: ";
		int size = 0;
		int min = 0;
		int max = 0;

		cin >> size;

		cout << "Enter maximum error: ";
		double error = 0;
		cin >> error;
		system("cls");

		matrix matr(size, size);
		matr.fillWithRandomForGradient(rand);
		mathVector vec(size);
		vec.fillWithRandom(rand);

		int uiRun = 1;
		int act = 0;
		while (uiRun)
		{

			cout << "1)Print matrix." << endl;
			cout << "2)Print vector." << endl;
			cout << "3)Enter matrix." << endl;
			cout << "4)Enter vector." << endl;
			cout << "5)Random matrix." << endl;
			cout << "6)Random vector." << endl;
			cout << "7)Gradient." << endl;
			cout << "8)Gauss." << endl;
			cout << "9)Stop." << endl;
			cout << endl << "Enter action: ";
			cin >> act;

			switch (act)
			{
			case 1:
				matr.printMatrix();

				system("pause");
				system("cls");
				break;
			case 2:
				vec.printVector();
				system("pause");
				system("cls");
				break;
			case 3:
				matr.fillWithHands();
				system("cls");
				break;
			case 4:
				vec.fillWithHands();
				system("cls");
				break;
			case 5:
				cout << "Enter minimum: ";
				cin >> min;
				cout << "Enter maximum: ";
				cin >> max;
				matr.fillWithRandomRangeForGradient(rand, min, max);
				system("cls");
				break;
			case 6:
				cout << "Enter minimum: ";
				cin >> min;
				cout << "Enter maximum: ";
				cin >> max;
				vec.fillWithRandomRange(rand, min, max);
				system("cls");
				break;
			case 7:
			{
				cout << "-------Gradient result-------" << endl;
				mathVector resultGradient = gradientMethod(matr, vec, error);
				cout << "Result gradient ";
				resultGradient.printVector();
				cout << "-----------------------------" << endl;
				system("pause");
				system("cls");
				break;
			}
			case 8:
			{
				cout << "-------Gauss result-------" << endl;
				mathVector resultGauss = gaussMethod(matr, vec, error);
				resultGauss.printVector();
				cout << "--------------------------" << endl;
				system("pause");
				system("cls");
				break;
			}
			case 9:
				uiRun = 0;
				system("cls");
				break;
			default:
				system("pause");
				system("cls");
				break;
			}
		}
		cout << "Run the program again? (0 - No/Another num - Yes)";
		cin >> progRun;
	}

}
