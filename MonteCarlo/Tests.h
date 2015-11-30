#pragma once
#ifndef TESTS_H
#define TESTS_H

#include <fstream>
#include <iostream>

using namespace std;

typedef double(*func1d) (double);
char csvSep = ',';

void test_1d(func1d f, double xMin, double xMax, unsigned nbPoints, char* fileName) {

	ofstream file;
	file.open(fileName);
	double step = (xMax - xMin) / (double)nbPoints;
	for (unsigned i = 0; i <= nbPoints; ++i)
		file << (double)i*step + xMin << csvSep << f((double)i*step + xMin) << endl;
	file.close();

}

template <typename Gen>
void progression_monte_carlo(long long n, Gen G, char fileName[] = "out.csv")
{
	ofstream file;
	file.open(fileName);
	double x = 0;
	for (long long j = 0; j < n; j++) {
		x += G();
		file << j+1 << csvSep << (x / (double)(j+1)) << endl;
		if (j%(long long)1e4 == 0)
			cout << ((double)j / (double)(n - 1)) * 100 << '%' << endl;
	}
	//file.close;
}

#endif \\ TESTS_H
