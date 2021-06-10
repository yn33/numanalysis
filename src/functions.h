#include "arctanx.cpp"
#include "lagrange.cpp"
#include "quadhermite.cpp"
#include "montecarlo.cpp"
#include "heun.cpp"
#include "fd.cpp"

double pi();

double arctanx(double x);

std::vector<double> lagweights(std::vector<double> xk);
std::vector<double> specialsum(std::vector<double> x, std::vector<double> z, std::vector<double> t);
std::vector<double> lagpolint(std::vector<double> t, std::vector<double> xk, std::vector<double> yk);
std::vector<double> quadhermite(double a, double b, double z, double y1, double y2, double s1, double s2);

double monteCarloBall(double n);

std::vector<double> form_t(int a, int b, double n);
std::vector<double> heun(std::vector<double> t, double (*fd)(double, double), double initial);

double fd_initial();
double fd(double t, double y);
double fd_solution(double t, double epsilon);

