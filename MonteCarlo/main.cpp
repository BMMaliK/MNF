#include "Tests.h"
#include "BS.h"
#include "monte-carlo.hpp"
#include "var_alea.hpp"

#include <iostream>
#include <fstream>

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

double logNormalize(double g) {
	double S0 = 100.;
	double r = .01;
	double sigma = .2;
	double T = 1.;

	return S0*exp((r - pow(sigma, 2) / 2) + sigma*g);
}

int main() {

	//test_1d(&callBS_Spot, 0., 200., 20000, "Call_BS.csv");

	// Prix Black Scholes
	cout << "Prix Black Scholes : " << callBS(100., 100., .2, 1., .01) << endl;

	// Convergence Monte Carlo
	init_alea();

	gaussian G;
	//progression_monte_carlo(2e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))), "Convergence_Monte_Carlo.csv");

	vector<double> result1e5, result1e6;
	ofstream file;
	file.open("Distribution_Monte_Carlo.csv");
	for (int i = 0; i < 100; ++i) {
		result1e5 = monte_carlo((long long)1e5, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))));
		result1e6 = monte_carlo((long long)1e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))));
		file << result1e5[0] << csvSep << result1e6[0] << endl;
		cout << "i = " << i << endl;
	}
	file.close();

	return 0;
}