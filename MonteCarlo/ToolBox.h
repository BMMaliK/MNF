#pragma once
#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <string>
#include "RandomGenerator.hpp"

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
	return x / (X.size() - 1);
};

double max(processus<double>::result_type X) {
	double x = 0;
	processus<double>::cst_iter i;
	for (i = ++X.begin(); i != X.end(); ++i)
		x = (x > (i->second)) ? x : i->second;
	return x;
};

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Pricing a European vanilla call option with a Monte Carlo method
double callMc(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, RandomGenerator g, string file_path = "out.csv") {
	ofstream fichier_gen(file_path, ios::out | ios::trunc);
	std::cout.precision(8);
	fichier_gen << "Nb_iterations" << "," << "Price" << endl;
	double spot = 0.0;
	double payoff_sum = 0.0;
	double spot_adjusted = S * exp(T*(r - 0.5*v*v));
	double spot_mc = 0.0;

	for (int i = 1; i<num_sims + 1; i++) {
		double gauss_bm = g.SuiteAleatoire();
		spot = spot_adjusted * exp(sqrt(v*v*T)*gauss_bm);
		payoff_sum += MAX(spot - K, 0.0);
		spot_mc = (payoff_sum / static_cast<double>(i)) * exp(-r*T);
		fichier_gen << i << "," << spot_mc << endl;
	}

	return spot_mc;
}

void WriteRandomToCsv(string file_path = "out.csv", int nb_points = 20000) {
	ofstream fichier_gen(file_path, ios::out | ios::trunc);


	fichier_gen << "Uniform" << "," << "Normal" << "," << "SQRT" << "," << "Halton" << endl;

	RandomGenerator Standard = RandomGenerator("Standard");
	RandomGenerator SQRT = RandomGenerator("SQRT");
	RandomGenerator Halton = RandomGenerator("Halton");

	for (int i = 0; i<nb_points; i++)
	{
		fichier_gen << Standard.uniform() << "," << Standard.SuiteAleatoire() << "," << SQRT.SuiteAleatoire() << "," << Halton.SuiteAleatoire() << endl;
	}

	fichier_gen.close();
}

#endif // TOOLBOX_H
