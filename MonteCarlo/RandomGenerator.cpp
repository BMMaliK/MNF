//
//  RandomGenerator.cpp
//  mnf_e
//

#include "RandomGenerator.hpp"
#include "math.h"
#define Pi 3.14159265358979323846
// We use class here of static contante or pointeurs because it is more elegant

RandomGenerator :: RandomGenerator()
{
    a = 16807;
    b = 0;
    M = pow(2.,31.) - 1;
    x=17;
    n=1;
    m=1;
    p=3;
    q=7;
    d=3;
    l=17;
}

RandomGenerator :: RandomGenerator(std::string NomMethod)
{
    this->NomMethod = NomMethod;
    a = 16807;
    b = 0;
    M = pow(2.,31.) - 1;
    x=17;
    n=1;
    m=1;
    p=5;
    q=19;
    d=5;
    l=7;
}
RandomGenerator :: ~RandomGenerator()
{
}

double RandomGenerator :: SuiteAleatoire()
{
    double s=1;
    if (NomMethod == "Standard")
        s = gaussian();
    else if (NomMethod == "Halton")
        s = Halton(d, l);
    else if (NomMethod == "SQRT")
        s = SQRT(p, q);
    return s;
}

double RandomGenerator :: uniform()
{
    x=(a*x+b)%M;
    return (double)x/M;
}
double RandomGenerator::gaussian()
{
    double u1;
    double u2;
    u1 = uniform();
    u2 = uniform();
    return sqrt((-2)*log(u1))*cos(2*Pi*u2);
}


double RandomGenerator :: SQRT(int p, int q)
{
    double U1;
    double U2;
    U1= n*sqrt((double) p) - floor(n*sqrt((double) p));
    U2= n*sqrt((double) q) - floor(n*sqrt((double) q));
    n+=1;
    //return U1;
    return sqrt(-2*log(U1))*cos(2*Pi*U2);
}


//HALTON
// Trouver l'activitÈ suivante
int RandomGenerator :: lemmeHalton(int p, int m)
{
    int i = 0;
    int exposantMax = 0;
    while ((m - pow((double)p,i)) >= 0 )
    {
        i++;
        exposantMax = i-1;
    }
    return exposantMax;
}
// Decomposition
double RandomGenerator :: halton(int p, int m)
{
    int exposant = lemmeHalton(p,m);
    int valo = m;
    double u = 0;
    while (valo!=0)
    {
        int coeff = valo / (int) pow((double)p,exposant);
        u += coeff / pow((double)p, exposant+1);
        valo = fmod (valo,pow((double)p,exposant));
        exposant = lemmeHalton(p, valo);
    }
    return u;
}
double RandomGenerator :: Halton(int d,int l)
{
    double U1;
    double U2;
    U1 = halton(d,m);
    U2 = halton(l,m);
    m += 1;
    //return U1;
    return sqrt(-2*log(U1))*cos(2*Pi*U2);;
}







