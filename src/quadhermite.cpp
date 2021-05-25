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