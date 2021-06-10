#include <math.h>

double fd_initial() {
    return 1.0/3.0;
}

double fd(double t, double y) {

    //this forms a differential equation
    double fdt = -150*y + 49 - 150*t;
    return fdt;
}

double fd_solution(double t, double epsilon) {
    
    //solution to the differential equation fd, for testing Heun/Euler
    double ft = -t + 1.0/3.0 + epsilon*exp(-150*t);
    return ft;
}