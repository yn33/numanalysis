#include <vector>
#include <exception>
#include <iostream>

// Forms the t vector from the given interval and n = 1/h
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

// Implements the Heun method given the vector t from form_t, the derivative function fd(t, y), and the initial value
// Returns a vector where every index corresponds to the index of the t value in the original vector t
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