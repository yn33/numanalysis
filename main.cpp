#include "src/errors.h"
#include <iostream>
#include <limits>

int main() {

    // Setting highest and lowest n, by default n increases logarithmically
    int nlimit = 10000000;
    int nstart = 100;

    // Ellipsoid dimensions for testing Monte Carlo
    int ellipsoid_a = 1;
    int ellipsoid_b = 2;
    int ellipsoid_c = 3;

    // Parameters for testing Heun's method
    int heun_a = 0;
    int heun_b = 1;
    int heun_epsilon = 0.01;

    // Setting the precision of values to show full double when printed
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    // Go through n logarithmically and run all error functions
    for(double n = nstart; n <= nlimit; n = n*10) {

        std::cout << "Current n:\n" << n << "\n";
        arctanError(n);
        lagrangeError(n);
        hermiteError(n);
        monteCarloError(n, ellipsoid_a, ellipsoid_b, ellipsoid_c);
        heunError(heun_a, heun_b, heun_epsilon, n);
        std::cout << "\n";
        
    }

    return 0;
}