#include <vector>
#include <iostream>
#include "functions.h"

// Error function for the arctan(x) approximation
// Divides pi/2 into n parts and evaluates these points

void arctanError(double n) {
    std::cout << "Arctan approximation\n";
    // Arctan approximation vs c++ atan()
    double maxerror = 0;
    double avgerror = 0;
    double pivalue = pi();
    for(int i = 0; i <= n; i++) {        
        double x = (pivalue/2)*(i/n);
        double cpptan = atan(x);
        double approxtan = arctanx(x);
        double error =  abs(cpptan - approxtan);
        avgerror = avgerror + error;
        if(error > maxerror) {
            maxerror = error;
        }
    }
    avgerror = avgerror/(n + 1);
    std::cout << "Maximum error:\n" << maxerror << "\n";
    std::cout << "Average error:\n" << avgerror << "\n";
}

// Error function for Lagrange interpolation
void lagrangeError(double n) {

    double pivalue = pi();
    std::cout << "Lagrange interpolation\n";
    // Lagrange interpolation for sqrt(abs(x)) in [-1, 1]
    // Uniform points
    std::vector<double> u{ -1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1 };
    std::vector<double> yu(9);
    for(int i = 0; i < 9; i++) {
        yu[i] = sqrt(abs(u[i]));
    }
    // t points
    std::vector<double> t(n + 1);
    for(int i = 0; i <= n; i++) {
        t[i] = -1 + 2*(i/n);
    }
    // Chebyshev points
    std::vector<double>yc(101);
    std::vector<double>xc(101);
    for(int i = 0; i <= 100; i++) {
        double xj = -cos((i/100.0)*pivalue);
        xc[i] = xj;
        yc[i] = sqrt(abs(xj));
    }
    std::vector<double> up = lagpolint(t, u, yu);
    std::vector<double> cp = lagpolint(t, xc, yc);
    double umaxerror = 0;
    double uavgerror = 0;
    double cmaxerror = 0;
    double cavgerror = 0;
    for(int i = 0; i <= n; i++) {
        double accurate = sqrt(abs(t[i]));
        double uest = up[i];
        double cest = cp[i];
        double uerror = abs(accurate - uest);
        double cerror = abs(accurate - cest);
        // The algorithm gives NaN for t = x, where error is 0
        if(isnan(uerror)) {
            uerror = 0;
        } 
        if(isnan(cerror)) {
            cerror = 0;
        }
        if(uerror > umaxerror) {
            umaxerror = uerror;
        }
        if(cerror > cmaxerror) {
            cmaxerror = cerror;
        }
        uavgerror = uavgerror + uerror;
        cavgerror = cavgerror + cerror;
    }
    uavgerror = uavgerror/(n + 1);
    cavgerror = cavgerror/(n + 1);
    std::cout << "Uniform points maximum error:\n" << umaxerror << "\n" << "Uniform points average error:\n" << uavgerror << "\n";
    std::cout << "Chebyshev points maximum error:\n" << cmaxerror << "\n" << "Chebyshev points average error:\n" << cavgerror << "\n";

}

// Error function for Hermite interpolation
void hermiteError(double n) {

    std::cout << "Quadratic hermite interpolation\n";
    // Quadratic hermite interpolation for x^4 on [0, 2]

    double a = 0;
    double b = 2;
    double z = 1;
    double y1 = 0;
    double y2 = 16;
    double s1 = 0;
    double s2 = 32;

    // Coefficients
    std::vector<double>c = quadhermite(a, b, z, y1, y2, s1, s2);
    std::vector<double>t_1(n + 1);
    std::vector<double>t_2(n + 1);

    // Equally distributing the vectors t on both sides of z
    for(int i = 0; i <= n; i++) {
        t_1[i] = (z - a)*(i/n);
        t_2[i] = (b - z)*(i/n);
    }
    
    // The interpolation
    std::vector<double>p = bernstein(n, t_1, t_2, c);

    double maxerror = 0;
    double avgerror = 0;

    for(int i = 0; i <= n; i++) {
        double currx = t_1[i];
        double accurate = pow(currx, 4);
        double interp = p[i];
        double error = abs(accurate - interp);

        if(error > maxerror) {
            maxerror = error;
        }
        avgerror = avgerror + error;
        currx = (z - a) + t_2[i];
        accurate = pow(currx, 4);
        interp = p[i + n + 1];
        error = abs(accurate - interp);
        if(error > maxerror) {
            maxerror = error;
        }
        avgerror = avgerror + error;
    }
    avgerror = avgerror/(2*(n + 1));

    std::cout << "Maximum error:\n" << maxerror << "\n";
    std::cout << "Average error:\n" << avgerror << "\n";
}

// Error function for Monte Carlo volume of ellipsoid
void monteCarloError(double n, int a, int b, int c) {

    std::cout << "Monte Carlo ellipsoid volume and pi approximations\n";

    double pivalue = pi();

    // Volumes of the hypercube and ellipsoid
    double cubeV = 2*a*2*b*2*c;
    double ellipsoidV = (4.0/3.0)*pivalue*a*b*c;
    
    // Approximation for the ratio of unit ball and unit cube volumes, pi/6
    double ballnk = monteCarloBall(n);

    // Ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 volume approximation
    double ellipsoidVapprox = cubeV*ballnk;
    double ellipsoidVerror = abs(ellipsoidV - ellipsoidVapprox);

    // Pi approximation
    double piapprox = ballnk*6;
    double pierror = abs(pivalue - piapprox);

    std::cout << "Ellipsoid volume error:\n" << ellipsoidVerror << "\n";
    std::cout << "Pi approximation:\n" << piapprox << "\n" << "Pi error:\n" << pierror << "\n";

}

// Error function for the Heun method
void heunError(int a, int b, double epsilon, double n) {

    std::cout << "Heun method\n";
    // Heun method with and without epsilon, the differential equation is set in fd.cpp

    double initial = fd_initial();

    // Set up t between a and b, with h = 1/n
    std::vector<double> t = form_t(a, b, n);
    double size = t.size();

    std::vector<double> u = heun(t, fd, initial);
    std::vector<double> u_epsilon = heun(t, fd, initial + epsilon);

    double uavgerror = 0;
    double umaxerror = 0;
    double uavgerror_epsilon = 0;
    double umaxerror_epsilon = 0;
    
    for(int i = 0; i < size; i++) {

        double o = fd_solution(t[i], 0);
        double o_epsilon = fd_solution(t[i], epsilon);
        double error = abs(o - u[i]);
        double error_epsilon = abs(o_epsilon - u_epsilon[i]);

        uavgerror = uavgerror + error;
        uavgerror_epsilon = uavgerror_epsilon + error_epsilon;
        if(error > umaxerror) {
            umaxerror = error;
        }
        if(error_epsilon > umaxerror_epsilon) {
            umaxerror_epsilon = error_epsilon;
        }
    }

    uavgerror = uavgerror/size;
    uavgerror_epsilon = uavgerror_epsilon/size;

    std::cout << "Maximum error without epsilon:\n" << umaxerror << "\n";
    std::cout << "Average error without epsilon:\n" << uavgerror << "\n";
    std::cout << "Maximum error with epsilon:\n" << umaxerror_epsilon << "\n";
    std::cout << "Average error with epsilon:\n" << uavgerror_epsilon << "\n";

}