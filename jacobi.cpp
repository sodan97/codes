#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/Core"
#include <fstream>
#include <chrono>

using namespace std;



Eigen::VectorXd jacobi(Eigen::MatrixXd& A, Eigen::VectorXd& b, int max_iterations, double tolerance) {
    int n = A.rows();
    Eigen::VectorXd U(n);
    U.setZero();

    for (int k = 0; k < max_iterations; k++) {
        Eigen::VectorXd U_new(n);
        U_new.setZero();

        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sigma += A(i, j) * U(j);
                }
            }
            U_new(i) = (b(i) - sigma) / A(i, i);
        }
        Eigen::VectorXd residual = b - A * U_new;
        double r_norm = residual.norm();
        U = U_new;

        if (r_norm < tolerance) {
            std::cout << "Convergence atteinte après " << k + 1 << " itérations." << std::endl;
            break;
        }
    }

    return U;
}


int main(){

    Eigen::VectorXd x;
    Eigen::MatrixXd A(2, 2);
    A <<  3, 2,
          2, 6;  
    Eigen::VectorXd b(2);
    b << 2, -8;
    int max_iterations = 1000;
    double tolerance = 1e-10;
    x = jacobi(A, b, max_iterations, tolerance);  
    cout << "A: \n" << A<< endl;
    cout << "b: \n" << b << endl;
    cout << "x: \n" << x << endl;

}





