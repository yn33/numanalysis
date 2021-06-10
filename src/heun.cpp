#include <vector>
#include <exception>
#include <iostream>

std::vector<double> form_t(int a, int b, double n) {

    double h = 1/n;
    int size = (b - a)*n + 1;
    std::vector<double> t(size);

    t[0] = a;

    for(int i = 1; i < size; i++) {
        t[i] = t[i - 1] + h;
    }

    int roundedb = round(t[size - 1]);
    if(roundedb != b) {
        std::cout << "Error: b was " << roundedb << " but it's supposed to be " << b << "\n";
        throw std::invalid_argument("Error when forming t at heun.cpp");
    }

    return t;
}

std::vector<double> heun(std::vector<double> t, double (*fd)(double, double), double initial) {

    double h = t[1] - t[0];
    int size = t.size();
    std::vector<double> u(size);
    double y = initial;
    u[0] = y;

    for(int i = 1; i < size; i++) {
        double prediction = y + h*fd(t[i - 1], y);
        y = y + (h/2)*(fd(t[i - 1], y), fd(t[i], prediction));
        u[i] = y;
    }

    return u;
}