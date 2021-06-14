#include <vector>
#include <exception>

// Computers the weights wk of the barycentric form Lagrange interpolation given the vector xk
std::vector<double> lagweights(std::vector<double> xk) {
    int n = xk.size();
    std::vector<double> wk(n);
    for(int i = 0; i < n; i++) {
        double currwk = 1;
        double currxk = xk[i];
        for(int j = 0; j < n; j++) {
            double currxi = xk[j];
            // Don't include i = j term
            if(i != j) {
                currwk = currwk*(1/(currxk - currxi));
            }
        }
        wk[i] = currwk;
    } 
    return wk;
}

// Computes the sum used in the barycentric form, where x is the vector containing the points where the Lagrange polynomial was constructed
// The vector z contains the weights
// The vector t contains the points where the interpolation is evaluated, and the result is a vector containing the sum for every given t
std::vector<double> specialsum(std::vector<double> x, std::vector<double> z, std::vector<double> t) {
    int n = x.size();
    int s = t.size();
    std::vector<double> res(s);
    if(n != z.size()) {
        throw std::invalid_argument("Wrong size input at specialsum(x, z, t)");
    }
    for(int i = 0; i < s; i++) {
        double currti = t[i];
        double sum = 0;
        for(int j = 0; j < n; j++) {
            double currxj = x[j];
            double currzj = z[j];
            sum = sum + currzj/(currti - currxj);
        }
        res[i] = sum;
    }
    return res;
}

// Forms the barycentric form Lagrange interpolation on the points of vector t
// Vector xk and yk are the points which the Lagrange polynomial is constructed from
std::vector<double> lagpolint(std::vector<double> t, std::vector<double> xk, std::vector<double> yk) {
    int n = xk.size();
    int s = t.size();
    std::vector<double> wk = lagweights(xk);
    if(n != yk.size()) {
        throw std::invalid_argument("Wrong size input at lagpolint(t, xk, yk)");
    }
    std::vector<double> res(s);
    std::vector<double> zk(n);
    // Forming the weights of the function values y
    for(int i = 0; i < n; i++) {
        zk[i] = wk[i]*yk[i];
    }
    std::vector<double> sum_1 = specialsum(xk, zk, t);
    std::vector<double> sum_2 = specialsum(xk, wk, t);
    // The final barycentric form result is the termwise division of these sums
    for(int i = 0; i < s; i++) {
        res[i] = sum_1[i]/sum_2[i];
    }
    return res;
}