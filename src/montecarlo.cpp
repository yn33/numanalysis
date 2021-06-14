#include <stdlib.h>
#include <math.h>
#include <iostream>

// This function forms three random numbers between 0 and 1, which are absolute values of coordinates inside an unit cube
// If the norm of these coordinates is less or equal to 1, they are contained within a unit ball
// This is repeated n times and used to approximate the volume  ratio between a unit ball and unit cube
double monteCarloBall(double n) {

    double nk = 0;
    for(int i = 0; i < n; i++) {
        double randmax = RAND_MAX;
        double randx = rand()/randmax;
        double randy = rand()/randmax;
        double randz = rand()/randmax;
        double radius = sqrt((pow(randx, 2) + pow(randy, 2) + pow(randz, 2)));

        if(radius <= 1.0) {
            nk++;
        }
    }
    return nk/n;
}