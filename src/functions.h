#include "arctanx.cpp"
#include "lagrange.cpp"
#include "quadhermite.cpp"
#include "montecarlo.cpp"

double pi();
double arctanx(double x);

std::vector<double> lagweights(std::vector<double> xk);
std::vector<double> specialsum(std::vector<double> x, std::vector<double> z, std::vector<double> t);
std::vector<double> lagpolint(std::vector<double> t, std::vector<double> xk, std::vector<double> yk);
std::vector<double> quadhermite(double a, double b, double z, double y1, double y2, double s1, double s2);
double monteCarloBall(double n);