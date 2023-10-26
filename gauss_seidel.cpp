#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/Core"
#include <fstream>
#include <chrono>

using namespace std;


Eigen::VectorXd gauss_seidel(Eigen::MatrixXd& A, Eigen::VectorXd& b, int max_iterations, double tolerance) {
    int n = A.rows();
    Eigen::VectorXd x(n), x_new(n);
    x.setZero();
    x_new.setZero();
    double error = tolerance + 1;
    int iter = 0; // Compteur d'itérations

    while (error > tolerance && iter < max_iterations) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                sum1 += A(i, j) * x_new(j);
            }

            for (int j = i + 1; j < n; j++) {
                sum2 += A(i, j) * x(j);
            }

            x_new(i) = (b(i) - sum1 - sum2) / A(i, i);
        }

        error = (x_new - x).norm(); // Calcul de l'erreur entre deux itérations
        x = x_new;

        if (error < tolerance) {
            cout << "Convergence atteinte après " << iter + 1 << " itérations." << endl;
            break;
        }
        iter++; // Incrémenter le compteur d'itérations
    }

    return x;
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
    x = gauss_seidel(A, b, max_iterations, tolerance);  
    cout << "A: \n" << A<< endl;
    cout << "b: \n" << b << endl;
    cout << "x: \n" << x << endl;

}