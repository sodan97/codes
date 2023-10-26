#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/Core"
#include <fstream>
#include <chrono>

using namespace std;


Eigen::VectorXd conjugate_gradient(Eigen::MatrixXd& A, Eigen::VectorXd& b, int max_iterations, double tolerance) {
    int n = A.rows();
    Eigen::VectorXd U(n), U_new(n);
    U.setZero(), U_new.setZero();
    double alpha;
    double error = tolerance + 1; 
    int iter = 0; // Compteur d'itérations
    while (error > tolerance && iter < max_iterations) {
        Eigen::VectorXd r = b - A * U; // Calcul du résidu  correspondant à la direction
        alpha = r.dot(r) / (r.dot(A * r)); // Calcul du pas (alpha) par la méthode de la descente de gradient
        U_new = U + alpha * r; // Mise à jour de la solution

        error = (U_new - U).norm();  // Calcul de l'erreur entre deux itérations
        U = U_new;

        ////////////////////    Conditions aux bords      //////////////////////
        // U_new(0) = 0;

        if (error < tolerance) {
            std::cout << "Convergence atteinte après " << iter + 1 << " itérations." << std::endl;
            break;
        }

        iter++; // Incrémenter le compteur d'itérations
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
    x = conjugate_gradient(A, b, max_iterations, tolerance);  
    cout << "A: \n" << A<< endl;
    cout << "b: \n" << b << endl;
    cout << "x: \n" << x << endl;

}