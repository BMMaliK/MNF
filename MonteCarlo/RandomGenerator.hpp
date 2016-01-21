//
//  RandomGenerator.hpp
//  mnf_e
//


#ifndef RandomGenerator_hpp
#define RandomGenerator_hpp
#include <string>
#include <stdio.h>

class RandomGenerator {
public:
    RandomGenerator ();
    RandomGenerator (std::string NomMethod);
    ~RandomGenerator();
    double SuiteAleatoire();
//    double boxmuller(double U1, double U2);
    double uniform();
    double gaussian();
    double SQRT(int p, int q);//SQRT

    //Halton
    double Halton(int d, int l);
    int lemmeHalton(int p, int m);
    double halton(int p, int m);

private:
    long long a, b, M;
public:
    std::string NomMethod ;
    long long x;
    int n,p,q;//Generateur SQRT
    int m,d,l;//Generateur Halton
};
#endif /* RandomGenerator_hpp */
