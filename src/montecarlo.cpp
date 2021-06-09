#include <stdlib.h>
#include <math.h>
#include <iostream>

double monteCarloBall(double n) {

    double nk = 0;
    for(int i = 0; i < n; i++) {
        double randmax = RAND_MAX;
        double randx = rand()/randmax;
        double randy = rand()/randmax;
        double randz = rand()/randmax;
        double radius = sqrt((pow(randx, 2) + pow(randy, 2) + pow(randz, 2)));
        double contained = 1.0 - radius;
        if(contained >= 0.0) {
            nk++;
        }
    }
    return nk/n;
}