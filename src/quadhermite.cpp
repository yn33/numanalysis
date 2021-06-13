#include <vector>

std::vector<double> quadhermite(double a, double b, double z, double y1, double y2, double s1, double s2) {
    double h1 = z - a;
    double h2 = b - z;

    double c10 = y1;                  
    double c22 = y2;                  
    double c11 = (h1/2)*s1 + y1;     
    double c21 = -(h2/2)*s2 + y2;     
    
    
    double c12 = (h1*c21 + h2*c11)/(h1+h2);
    double c20 = c12;

    std::vector<double>c { c10, c11, c12, c20, c21, c22 };
    return c;
}

std::vector<double> bernstein(double n, std::vector<double> t_1, std::vector<double> t_2, std::vector<double> c) {

    std::vector<double>B_10(n + 1);
    std::vector<double>B_11(n + 1);
    std::vector<double>B_12(n + 1);
    std::vector<double>B_20(n + 1);
    std::vector<double>B_21(n + 1);
    std::vector<double>B_22(n + 1);

    for(int i = 0; i <= n; i++) {
        B_10[i] = pow((1 - t_1[i]), 2);
        B_20[i] = pow((1 - t_2[i]), 2);
        B_11[i] = 2*t_1[i]*(1 - t_1[i]);
        B_21[i] = 2*t_2[i]*(1 - t_2[i]);
        B_12[i] = pow(t_1[i], 2);
        B_22[i] = pow(t_2[i], 2);
    }
    
    std::vector<double>p(2*(n + 1));
    for(int i = 0; i <= n; i++) {
        p[i] = c[0]*B_10[i] + c[1]*B_11[i] + c[2]*B_12[i];
        p[i + n + 1] = c[3]*B_20[i] + c[4]*B_21[i] + c[5]*B_22[i];
    }

    return p;
}