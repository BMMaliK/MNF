#include "BinTree.h"

double BinTree::payoff()
{
	return (value >= K) ? (value - K) : 0;
}

double BinTree::delta(double relative)
{
	return (getShiftedValueClone(relative)->recombine() - recombine())/(value*relative);
}

double BinTree::gamma(double relative)
{
	return (delta(relative)-delta(-relative))/(value*relative);
}

double BinTree::vega(double relative)
{
	return (getShiftedVolClone(relative)->recombine() - recombine())/(vol*relative);
}

double BinTree::theta(double relative)
{
	return (getShiftedMaturityClone(relative)->recombine() - recombine()) / (T*relative);
}

double BinTreePut::payoff()
{
	return (value < K) ? -(value - K) : 0;
}
