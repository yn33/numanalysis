#include "functions.h"
#include <vector>
#include <iostream>

int main() {
    //arctan approximation vs c++ atan()
    double maxerror = 0;
    double avgerror = 0;
    double n = 10000;
    double pivalue = pi();
    for(int i = 0; i < n; i++) {        
        double x = (pivalue/2)*(i/n);
        double cpptan = atan(x);
        double approxtan = arctanx(x);
        double error =  abs(cpptan - approxtan);
        avgerror = avgerror + error;
        if(error > maxerror) {
            maxerror = error;
        }
    }
    avgerror = avgerror/n;
    std::cout << "Maximum error: " << "\n";
    std::cout << maxerror << "\n";
    std::cout << "Average error: " << "\n";
    std::cout << avgerror << "\n";
    //maxerror 1.1191e-13, pretty good approximation
    //avgerror 2.00266e-14
    return 0;
}