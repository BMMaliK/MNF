#pragma once
#ifndef BINTREE_H
#define BINTREE_H

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;

class BinTree
{
public:
	BinTree(double value, unsigned int depth, double vol, double K, double r = 0., double T = 1., BinTree* down = nullptr, bool build = true) : depth(depth), optionValue(0.), calculated(false),  down(down), value(value), T(T), step(T / (double)depth), r(r), vol(vol), K(K) {
		if(build) this->build();
	};
	BinTree* getUp() { return up; };
	double u() { return exp(vol*sqrt(step)); };
	double d() { return exp(-vol*sqrt(step)); };
	double p() { return (exp(r*step) - d()) / (u() - d()); };
	bool isLeaf() { return (down == nullptr); };
	virtual double recombine() {
		if (isLeaf())
			optionValue = payoff();
		else {
			if (!calculated) {
				optionValue = exp(-r*step)*(p()*up->recombine() + (1 - p())*down->recombine());
			}
		}
		calculated = true;
		return optionValue;
	};
	virtual double payoff();
	void clear() {	
		if (calculated) {
			calculated = false;
			optionValue = 0.;
			if (!isLeaf()) {
				up->clear();
				down->clear();
			}
		}
	};
	double delta(double relative = .01);
	double gamma(double relative = .01);
	double vega(double relative = .01);
	double theta(double relative = .01);

protected:
	double value; //S, Node Value
	double T; //Maturity
	double step;
	double r; //Interest rate
	double vol; //Volatility
	double K; //Option Strike
	unsigned int depth;

	bool calculated;
	double optionValue;

	BinTree* down;
	BinTree* up;

	virtual void build() {
		if (depth != 0) {
			if (down == nullptr) {
				down = new BinTree(value*d(), depth - 1, vol, K, r, T - step);
			}
			if (down == nullptr)
				up = new BinTree(value*u(), depth - 1, vol, K, r, T - step);
			else
				up = new BinTree(value*u(), depth - 1, vol, K, r, T - step, this->down->getUp());
		}
		else {
			this->down = nullptr;
			up = nullptr;
		}
	};
	virtual BinTree* clone() {
		return new BinTree(value,depth,vol,K,r,T,nullptr,false);
	};
	BinTree* getShiftedValueClone(double relative) {
		BinTree* cloned = clone();
		cloned->value = value*(1 + relative);
		cloned->build();
		return cloned;
	};
	BinTree* getShiftedVolClone(double relative) {
		BinTree* cloned = clone();
		cloned->vol = vol*(1 + relative);
		cloned->build();
		return cloned;
	};
	BinTree* getShiftedMaturityClone(double relative) {
		BinTree* cloned = clone();
		cloned->T = T*(1 + relative);
		cloned->build();
		return cloned;
	};
};

class BinTreePut : public BinTree {
public:
	BinTreePut(double value, unsigned int depth, double vol, double K, double r = 0., double T = 1., BinTree* down = nullptr, bool build = true) : BinTree(value, depth, vol, K, r, T, down, false) {
		if(build) this->build();
	};

protected:
	void build() {
		if (depth != 0) {
			if (down == nullptr) {
				this->down = new BinTreePut(value*d(), depth - 1, vol, K, r, T - step);
			}
			if (this->down == nullptr)
				up = new BinTreePut(value*u(), depth - 1, vol, K, r, T - step);
			else
				up = new BinTreePut(value*u(), depth - 1, vol, K, r, T - step, this->down->getUp());
		}
		else {
			this->down = nullptr;
			up = nullptr;
		}
	};
	virtual BinTree* clone() {
		return new BinTreePut(value, depth, vol, K, r, T, nullptr, false);
	};
	double payoff();
};

class AmericanBinTreePut : public BinTreePut {
public:
	AmericanBinTreePut(double value, unsigned int depth, double vol, double K, double r = 0., double T = 1., BinTree* down = nullptr, bool build = true) : BinTreePut(value, depth, vol, K, r, T, down, false) {
		if (build) this->build();
	};
	double recombine() {
		optionValue = BinTree::recombine();
		optionValue = (optionValue > payoff()) ? optionValue : payoff();
		return optionValue;
	};

protected:
	void build() {
		if (depth != 0) {
			if (down == nullptr) {
				this->down = new AmericanBinTreePut(value*d(), depth - 1, vol, K, r, T - step);
			}
			if (this->down == nullptr)
				up = new AmericanBinTreePut(value*u(), depth - 1, vol, K, r, T - step);
			else
				up = new AmericanBinTreePut(value*u(), depth - 1, vol, K, r, T - step, this->down->getUp());
		}
		else {
			this->down = nullptr;
			up = nullptr;
		}
	};
	virtual BinTree* clone() {
		return new AmericanBinTreePut(value, depth, vol, K, r, T, nullptr, false);
	};
};

class AmericanBinTree : public AmericanBinTreePut {
public:
	AmericanBinTree(double value, unsigned int depth, double vol, double K, double r = 0., double T = 1., BinTree* down = nullptr, bool build = true) : AmericanBinTreePut(value, depth, vol, K, r, T, down, false) {
		if (build) this->build();
	};

protected:
	void build() {
		if (depth != 0) {
			if (down == nullptr) {
				this->down = new AmericanBinTree(value*d(), depth - 1, vol, K, r, T - step);
			}
			if (this->down == nullptr)
				up = new AmericanBinTree(value*u(), depth - 1, vol, K, r, T - step);
			else
				up = new AmericanBinTree(value*u(), depth - 1, vol, K, r, T - step, this->down->getUp());
		}
		else {
			this->down = nullptr;
			up = nullptr;
		}
	};
	double payoff() { return BinTree::payoff(); };
	virtual BinTree* clone() {
		return new AmericanBinTree(value, depth, vol, K, r, T, nullptr, false);
	};
};


#endif //BINTREE_H

