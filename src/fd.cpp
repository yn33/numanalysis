#include <math.h>

// Initial value without error term
double fd_initial() {
    return 1.0/3.0;
}

// The derivative function which forms the differential equation
double fd(double t, double y) {

    double fdt = -150*y + 49 - 150*t;
    return fdt;
}

// Solution to the differential equation fd, for testing Heun/Euler
double fd_solution(double t, double epsilon) {
    
    double ft = -t + 1.0/3.0 + epsilon*exp(-150*t);
    return ft;
}