#include "src/errors.h"
#include <iostream>
#include <limits>

int main() {

    int nlimit = 10000000;
    int a = 1;
    int b = 2;
    int c = 3;

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    for(double n = 1; n <= nlimit; n = n*10) {

        std::cout << "Current n:\n" << n << "\n";
        arctanError(n);
        lagrangeError(n);
        hermiteError(n);
        monteCarloError(n, a, b, c);
    }

    return 0;
}