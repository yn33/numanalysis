#include <math.h>

double pi() {
    return 3.141592653589793;
}

double arctanx(double x) {
    double pi = 3.141592653589793;
    double s;
    double y;
    double a;
    double b;
    double c;
    double d;
    double u;
    double atanu;
    if (x >= 0.0 && x <= 1.7*pow(10, -9)) {
        s = x;
    }
    else if (x > 1.7*pow(10, -9) && x <= 2*pow(10, -2)) {
        s = x - pow(x, 3)/3 + pow(x, 5)/5 - pow(x, 7)/7;
    } else {
        if (x >= 0 && x <= 1) {
            y = x;
            a = 0;
            b = 1;
        }
        if (x > 1) {
            y = 1/x;
            a = pi/2;
            b = -1;  
        }
        c = pi/16;       
        if (y > (sqrt(2)-1) && y <= 1) {
            c = (3*pi)/16;
        }
        d = tan(c);
        u = (y - d)/(1 + d*y);
        atanu = u*((135135+171962.46*pow(u, 2)+52490.4832*pow(u, 4)+2218.1*pow(u, 6))/(135135+217007.46*pow(u, 2)+97799.3033*pow(u, 4)+10721.3745*pow(u, 6)));
        s = a + b*(c + atanu);
    }
    return s;
}