#include "Tests.h"
#include "BS.h"
#include "monte-carlo.hpp"
#include "var_alea.hpp"
#include "processus.h"

#include <iostream>
#include <fstream>

using namespace std;

double actualize(double intrinsic) {
	double r = .01;
	double T = 1.;

	return exp(-r*T)*intrinsic;
}

double payoff(double S_T) {
	//Call
	double K = 100.;
	return (S_T > K) ? (S_T - K) : 0;
}

double antitheticPayoff(pair<double, double> intrinsic) {
	return .5*(payoff(intrinsic.first) + payoff(intrinsic.second));
}

double logNormalize(double g) {
	double S0 = 100.;
	double r = .01;
	double sigma = .2;
	double T = 1.;

	return S0*exp((r - pow(sigma, 2) / 2) + sigma*g);
}

std::pair<double, double> antitheticLogNormalize(double g) {
	pair<double, double> temp(logNormalize(g), logNormalize(-g));
	return temp;
}

processus<double>::result_type logNormalizeProc(processus<double>::result_type X) {
	double S0 = 100.;
	double r = .01;
	double sigma = .2;
	double T = 1.;
	processus<double>::result_type value(X.size());

	processus<double>::cst_iter i;
	processus<double>::iter j;
	for (i = X.begin(), j = value.begin(); i != X.end(); ++i, ++j) {
		j->first = i->first;
		j->second = S0*exp((r - pow(sigma, 2) / 2)*i->first + sigma*i->second);
	}
	return value;
}

double average(processus<double>::result_type X) {
	double x = 0;
	processus<double>::cst_iter i;
	for (i = ++X.begin(); i != X.end(); ++i)
		x += i->second;
	return x / (X.size()-1);
};

double max(processus<double>::result_type X) {
	double x = 0;
	processus<double>::cst_iter i;
	for (i = ++X.begin(); i != X.end(); ++i)
		x = (x > (i->second)) ? x : i->second;
	return x;
};

int main() {

	//test_1d(&callBS_Spot, 0., 200., 20000, "Call_BS.csv");

	// Prix Black Scholes
	cout << "Prix Black Scholes : " << callBS(100., 100., .2, 1., .01) << endl;

	// Convergence Monte Carlo
	init_alea();

	//gaussian G;
	////progression_monte_carlo(2e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))), "Convergence_Monte_Carlo.csv");

	vector<double> result1e5, result1e6;
	//ofstream file;
	//file.open("Distribution_Monte_Carlo.csv");
	//for (int i = 0; i < 100; ++i) {
	//	result1e5 = monte_carlo((long long)1e5, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))));
	//	result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))));
	//	file << result1e5[0] << csvSep << result1e6[0] << endl;
	//	cout << "i = " << i << endl;
	//}
	//file.close();

	//Standard
	gaussian G;
	result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))));
	cout << "Standard : Prix = " << result1e6[0] << endl << "Ecart-Type = " << sqrt(result1e6[1]/1e6) << endl;
	result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(antitheticPayoff), compose(ptr_fun(antitheticLogNormalize), G))));
	cout << "Antithetique : Prix = " << result1e6[0] << endl << "Ecart-Type = " << sqrt(result1e6[1]/1e6) << endl;

	//Asian
	unsigned int N = 100;
	brownian B(N);
	result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(average),compose(ptr_fun(logNormalizeProc), B)))));
	cout << "Prix avec asianing sur N = " << N << " : " << result1e6[0] << endl;


	//LookBack
	result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(max), compose(ptr_fun(logNormalizeProc), B)))));
	cout << "Prix avec lookBack sur N = " << N << " : " << result1e6[0] << endl;

	return 0;
}