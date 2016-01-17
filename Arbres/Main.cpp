#include <iostream>

#include "BinTree.h"

using namespace std;

int main() {
	double S0 = 100.;
	unsigned int depth = 1000;
	double vol = .2;
	double K = 100.;
	double r = .01;
	double T = 1.;

	BinTree arbreCall(S0, depth, vol, K, r, T);
	BinTreePut arbrePut(S0, depth, vol, K, r, T);
	AmericanBinTreePut arbrePutAm(S0, depth, vol, K, r, T);
	AmericanBinTree arbreCallAm(S0, depth, vol, K, r, T);

	double prixCall = arbreCall.recombine();
	double prixPut = arbrePut.recombine();
	double prixPutAm = arbrePutAm.recombine();
	double prixCallAm = arbreCallAm.recombine();

	double deltaCall = arbreCall.theta();
	double deltaPut = arbrePut.theta();
	double deltaPutAm = arbrePutAm.theta();
	double deltaCallAm = arbreCallAm.theta();

	cout << "Prix Call = " << prixCall << endl;
	cout << "Prix Put = " << prixPut << endl;
	cout << "Prix Put Am = " << prixPutAm << endl;
	cout << "Prix Call Am = " << prixCallAm << endl;

	cout << "Delta Call = " << deltaCall << endl;
	cout << "Delta Put = " << deltaPut << endl;
	cout << "Delta Put Am = " << deltaPutAm << endl;
	cout << "Delta Call Am = " << deltaCallAm << endl;

	cout << prixCall - prixPut << endl;
	cout << 100.*(1 - exp(-.01*1.)) << endl;
	return 0;
}