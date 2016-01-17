#include "Tests.h"
#include "BS.h"
#include "monte-carlo.hpp"
#include "var_alea.hpp"
#include "processus.h"
#include "ToolBox.h"
#include "RandomGenerator.hpp"

#include <iostream>
#include <fstream>

using namespace std;

int main() {
	RandomGenerator Standard = RandomGenerator("Standard");
	RandomGenerator SQRT = RandomGenerator("SQRT");
	RandomGenerator Halton = RandomGenerator("Halton");
	int num_sims = 1e6;
	double S = 100.0;  // Option price
	double K = 100.0;  // Strike price
	double r = 0.01;   // Risk-free rate (1%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.0;
	//double call_price = callBS(S, K, v, T, r);
	//double call_price_mc = callMc(num_sims, S, K, r, v, T, Standard);
	//callMc(num_sims, S, K, r, v, T, SQRT, "mc_sqrt.csv");
	//callMc(num_sims, S, K, r, v, T, Halton, "mc_halton.csv");

	// Tracer le prix du Call
	//test_1d(&callBS_Spot, 0., 200., 20000, "Call_BS.csv");

	// Prix Black Scholes
	cout << "Prix Black Scholes : " << callBS(S, K, v, T, r) << endl;

	init_alea();

	// Convergence Monte Carlo
	//gaussian G;
	////progression_monte_carlo(2e6, compose(ptr_fun(actualize), compose(ptr_fun(payoff), compose(ptr_fun(logNormalize), G))), "Convergence_Monte_Carlo.csv");

	vector<double> result1e6;
	//Densité
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