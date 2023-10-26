#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/Core"
#include <fstream>
#include <chrono>

using namespace std;


//////////////////////////////////////////////////////////////////    Cas Poutre 1D       ////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                      //
//                           Chaque noeud a deux degrés de liberté (suivant x: traction ou compression; rotation autour de l'axe neutre)                                //
//                                                                                                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fonction pour calculer la matrice de rigidité locale pour la poutre 1D
Eigen::MatrixXd local_stiffness_matrix_beam1D(double x1, double x2,  double E, double I) {
    double L = abs(x2 - x1); 
    Eigen::MatrixXd Ke(4,4);
    Ke.setZero();
    Ke(0,0) = 12;       Ke(0,1) = 6*L;          Ke(0,2) = -12;          Ke(0,3) = 6*L;
    Ke(1,0) = 6*L;      Ke(1,1) = 4*L*L;        Ke(1,2) = -6*L;         Ke(1,3) = 2*L*L;
    Ke(2,0) = -12;      Ke(2, 1) = -6*L;        Ke(2,2) = 12;           Ke(2,3) = -6*L;
    Ke(3,0) = 6*L;      Ke(3,1) = 2*L*L;        Ke(3,2) = -6*L;         Ke(3,3) = 4*L*L; 
    Ke = (E*I/pow(L,3)) * Ke;
    return Ke;
}



Eigen::MatrixXd global_stiffness_matrix_beam1D(double length, double E, double I, int nbOfBeam) {
    int nbNodes = nbOfBeam + 1; 
    Eigen::MatrixXd K(2 * nbNodes , 2 * nbNodes);
    double delta_x = length/nbOfBeam; 
    K.setZero();

    double x1, x2;

    for (int ind_beam = 0; ind_beam < nbOfBeam; ind_beam++) {
        x1 = ind_beam * delta_x;
        x2 = (ind_beam + 1) * delta_x;
        Eigen::MatrixXd Ke = local_stiffness_matrix_beam1D(x1, x2, E, I);

        // Assemblage de la matrice locale dans la matrice globale
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                K(2* ind_beam + i, 2 * ind_beam + j) += Ke(i, j);
            }
        }
    }

    return K;
}

// Eigen::VectorXd global_load_vector(int nbOfBeam){
//     int nbNodes = nbOfBeam + 1; 
//     Eigen::VectorXd  F(2 * nbNodes);
//     F.setZero();
//     int indice = F.size() - 2;
//     F(indice) = 4000; //100000;                   // en Newton
//     return F;
// }


//////////////////////////////////////////////////////////////////    Cas Poutre 2D       ////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                      //
//                           Chaque noeud a trois degrés de liberté (suivant x: traction ou compression; suivant y: flexion; rotation autour de l'axe neutre)           //
//                                                                                                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fonction pour calculer la matrice de rigidité locale pour la poutre 2D
Eigen::MatrixXd local_stiffness_matrix_beam2D(double x1, double x2, double A, double C, double E, double I, double S) {
    double L = abs(x2 - x1); 
    Eigen::MatrixXd Ke(6,6);
    Ke.setZero();
    Ke(0,0) = A*C*C + 12*I*S*S/(L*L);  Ke(0,1) = (A-12*I/(L*L))*C*S;       Ke(0,2) = -6*I*S/L;   Ke(0,3) = -(A*C*C + 12*I*S*S/(L*L));  Ke(0,4) = -(A-12*I/(L*L))*C*S;    Ke(0,5) = -6*I*S/L; 

                                       Ke(1,1) = A*S*S + 12*I*C*C/(L*L);   Ke(1,2) = 6*I*C/L;    Ke(1,3) = -Ke(0,1);                   Ke(1,4) = -Ke(1,1);               Ke(1,5) = Ke(1,2);

                                                                           Ke(2,2) = 4*I;        Ke(2,3) = -Ke(0,2);                   Ke(2,4) = -Ke(1,2);               Ke(2,5) = 2*I;
                                                                                                
                                                                                                 Ke(3,3) = -Ke(0,3);                   Ke(3,4) = Ke(0,4);                Ke(3,5) = -Ke(0,5);
        
                                                                                                                                       Ke(4,4) = Ke(1,1);                Ke(4,5) = -Ke(1,2);
                                                        
                                                                                                                                                                         Ke(5,5) = Ke(2,2);

    for (int i=1; i<6; i++){
        for (int j=0; j<i; j++){
            Ke(i,j) = Ke(j,i);
        }
    }
    Ke = (E/L) * Ke;
    return Ke;
}



Eigen::MatrixXd global_stiffness_matrix_beam2D(double length, double A, double C, double E, double I, double S, int nbOfBeam) {
    int nbNodes = nbOfBeam + 1; 
    Eigen::MatrixXd K(3 * nbNodes ,3 * nbNodes);
    double delta_x = length/nbOfBeam; 
    K.setZero();

    double x1, x2;

    for (int ind_beam = 0; ind_beam < nbOfBeam; ind_beam++) {
        x1 = ind_beam * delta_x;
        x2 = (ind_beam + 1) * delta_x;
        Eigen::MatrixXd Ke = local_stiffness_matrix_beam2D(x1, x2, A, C, E, I, S);

        // Assemblage de la matrice locale dans la matrice globale
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                K(3* ind_beam + i, 3 * ind_beam + j) += Ke(i, j);
            }
        }
    }

    return K;
}



// Eigen::VectorXd global_load_vector(int nbOfBeam){
//     int nbNodes = nbOfBeam + 1; 
//     Eigen::VectorXd  F(3 * nbNodes);            // Force 
//     F.setZero();       

//     // Cas flexion             
//     // int indice = F.size() - 2;
//     // F(indice) = 4000; //100000;      

//     // Cas traction
//     int indice = F.size() - 3;        
//     F(indice) = 10000;              
//     return F;
// }



//////////////////////////////////////////////////////////////////    Cas Poutre 3D       ////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                      //
//                           Chaque noeud a six degrés de liberté: ux, theta_x                                                                                          //
//                                                                 uy, theta_y                                                                                          //
//                                                                 uz, theta_z                                                                                          //
//                                                                                                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Eigen::MatrixXd transformation_matrix(Eigen::MatrixXd coords, double length){

    double l = (coords(1,0) - coords(0,0))/length;                   //x2 = coords(1,0) et x1 = coords(0,0)
    double m = (coords(1,1) - coords(0,1))/length;                   //y2 = coords(1,1) et y1 = coords(0,1)
    double n = (coords(1,2) - coords(0,2))/length;                   //z2 = coords(1,2) et z1 = coords(0,2)
    double D = pow(l*l  + m*m, 0.5);
    Eigen::MatrixXd T_sequence(3, 3);
    Eigen::MatrixXd T(12, 12);
    T_sequence.setZero(); T.setZero();
    T_sequence <<  l,           m,          n,
                  -m/D,        l/D,         0,
                  -l*n/D,      -m*n/D,      D;

  // Copiez T_sequence dans la diagonale supérieure droite de T
    T.block(0, 0, 3, 3) = T_sequence;

   // Copiez T_sequence dans la diagonale inférieure droite de T
    T.block(3, 3, 3, 3) = T_sequence;

   // Copiez T_sequence dans la diagonale inférieure droite de T
    T.block(6, 6, 3, 3) = T_sequence;

   // Copiez T_sequence dans la diagonale inférieure droite de T
    T.block(9, 9, 3, 3) = T_sequence;

    return T;

}


// Fonction pour calculer la matrice de rigidité locale dans le cas de la poutre 3D
Eigen::MatrixXd local_stiffness_matrix_beam3D( double x1, double x2, Eigen::MatrixXd coords, double A, double E, double Iy, double Iz,double G, double J, double length) {
    double L = abs(x2 - x1);                                                                            //
    // double L = length;                                                                               /*/
    Eigen::MatrixXd Ke(12,12); 
    Eigen::MatrixXd T(12,12); 
    Ke.setZero(); T.setZero();

    Ke(0,0) = A*E/L; Ke(0,6) = -Ke(0,0) ; Ke(1,1) = 12*E*Iz/(pow(L,3));  Ke(1,5) = 6*E*Iz/(pow(L,2)); Ke(1,7) = -Ke(1,1); Ke(1,11) = Ke(1,5);

    Ke(2,2) = 12*E*Iz/(pow(L,3));  Ke(2,4) = -6*E*Iy/(pow(L,2)); Ke(2,8) = -Ke(2,2);  Ke(2,10) = Ke(2,4);  Ke(3,3) = G*J/L; Ke(3,9) = -Ke(3,3); 

    Ke(4,4) = 4*E*Iy/L;  Ke(4,8) = 6*E*Iy/(pow(L,2));    Ke(4,10) =2*E*Iy/L;
    
    Ke(5,5) = 4*E*Iz/L;  Ke(5,7) = -6*E*Iz/(pow(L,2));   Ke(5,11) = 2*E*Iz/L;

    Ke(6,6) = Ke(0,0); Ke(7,7) = Ke(1,1); Ke(7,11) = Ke(5,7);  Ke(8,8) = 12*E*Iy/(pow(L,3));   Ke(8,10) = Ke(4,8); 

    Ke(9,9) = Ke(3,3); Ke(10,10) = Ke(4,4);  Ke(11,11) = Ke(5,5);

    for (int i=1; i<12; i++){
        for (int j=0; j<i; j++){
            Ke(i,j) = Ke(j,i);
        }
    }

    T = transformation_matrix(coords, length);
    Ke = T.transpose() * Ke * T;

    return Ke;
}



Eigen::MatrixXd global_stiffness_matrix_beam_3D(Eigen::MatrixXd coords, double A, double E, double Iy, double Iz,double G, double J, double length, int nbOfBeam) {
    int nbNodes = nbOfBeam + 1; 
    Eigen::MatrixXd K(6 * nbNodes, 6 * nbNodes);
    double delta_x = length/nbOfBeam; 
    K.setZero();

    double x1, x2;                                                                                                         //

    for (int ind_beam = 0; ind_beam < nbOfBeam; ind_beam++) {
        x1 = ind_beam * delta_x;                                                                                        //
        x2 = (ind_beam + 1) * delta_x;                                                                                  //
        Eigen::MatrixXd Ke = local_stiffness_matrix_beam3D(x1, x2, coords, A, E, Iy, Iz, G, J, length);

        // Assemblage de la matrice locale dans la matrice globale
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                K(6* ind_beam + i, 6 * ind_beam + j) += Ke(i, j);
            }
        }
    }

    return K;
}







//////////////////////////////////////////////////////////////////    FORCE        /////////////////////////////////////////////////////////////////////////////////////


Eigen::VectorXd global_load_vector(int nbOfBeam){
    int nbNodes = nbOfBeam + 1; 
    Eigen::VectorXd  F(6 * nbNodes);            // Force 
    F.setZero();       

    // Cas flexion             
    int indice = F.size() - 4;
    F(indice) = 4000; //100000;      

    // // Cas traction
    // int indice = F.size() - 6;        
    // F(indice) = 21360; //10000;       


    return F;
}



///////////////////////////////////////////////////////////////// Methode du gradient conjugué: ////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd conjugate_gradient(Eigen::MatrixXd& A, Eigen::VectorXd& b, int max_iterations, double tolerance) {
    ofstream file("ResultatsFlexion.txt");
    int n = A.rows();
    Eigen::VectorXd U(n), U_new(n);
    U.setZero(), U_new.setZero();
    double alpha;
    double error = tolerance + 1; 
    int iter = 0;                                           // Compteur d'itérations
    while (error > tolerance && iter < max_iterations) {
        Eigen::VectorXd r = b - A * U;                      // Calcul du résidu  correspondant à la direction
        alpha = r.dot(r) / (r.dot(A * r));                  // Calcul du pas (alpha) par la méthode de la descente de gradient
        U_new = U + alpha * r;                              // Mise à jour de la solution

        ////////////////////    Conditions aux bords      //////////////////////
        U_new(0) = 0;
        U_new(1) = 0;
        U_new(2) = 0;
        U_new(3) = 0;
        U_new(4) = 0;
        U_new(5) = 0;

        error = (U_new - U).norm();  // Calcul de l'erreur entre deux itérations
        U = U_new;
        // file << iter << "        " << U(n-3) << endl;          // traction
        // file << iter << "        " << U(n-2) << endl;          // flexion 2D
        file << iter << "        " << U(n-4) << endl;          // flexion 3D
        if (error < tolerance) {
            std::cout << "Convergence atteinte après " << iter + 1 << " itérations." << std::endl;
            break;
        }
        
        iter++; // Incrémenter le compteur d'itérations
    }
    return U;
}





int main() {

    ofstream file2("exactSolution.txt");
    for (int i = 0; i <= 100000; i++){
        file2 << i << "   " << 1.48e-2 << endl;                // flexion
        //  file2 << i << "   " << 666e-6 << endl;             // traction
    }
    auto init = chrono::high_resolution_clock::now();

    /////////////////////////////////////////////     Cas flexion 3D    /////////////////////////////////////////////

    double length = 2; //120; //2.0; //5; // Longueur totale de la poutre (m)
    double E = 200e9;  // 72e9;  // Module de Young (en Gpa)
    double I = 3.6e-6; //200; //3.6e-6; // 1e-8/12;//3.6e-6;  //4e-5; // Moment d'inertie (m^4)
    int nbOfBeam = 5;   // Nombre d'éléments de poutre
    int max_iterations = 100000;
    double tolerance = 1e-7;
    double A = 200e-3 * 60e-3; //10;    // section en mètre carré

    Eigen::MatrixXd coords(2,3);
    // double Iy = 0; double Iz = 0; double G = 10000; double J = 50; 
    double Iy = 3.6e-6; double Iz = 3.6e-6; double G = 73e9; double J = 0;   // J: constante de torsion pour les sections non circulaire et moment d'inertie polaire pour les sections circulaires
    coords << 0,   0,   0,
              2,   0,   0;

    Eigen::MatrixXd K_global = global_stiffness_matrix_beam_3D(coords, A, E, Iy, Iz, G, J, length, nbOfBeam);
    
    Eigen::VectorXd  F = global_load_vector(nbOfBeam);                                                         
                  
    // Eigen::SparseMatrix<double> K_global_sparse = K_global.sparseView(); // Conversion en format sparse
    // cout << "K_global: \n" << K_global << endl;
    // cout << " K_global_sparse: \n" << K_global_sparse << endl; 

    Eigen::VectorXd U = conjugate_gradient(K_global, F, max_iterations, tolerance);                              
    // Eigen::VectorXd U = conjugate_gradient(K_global_sparse, F, max_iterations, tolerance);


    // Afficher la solution U
    std::cout << "Solution U :\n" << U << std::endl;                                                     

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - init).count();
    cout << "Temps d'éxécution (en secondes): " << duration << endl;

    return 0;
}











































///////////////////////////////////////////////////////////////////       Poutre 1D | 2D      ///////////////////////////////////////////////////////////////////////////////


// int main() {

//     ofstream file2("exactSolution.txt");
//     for (int i=0; i<=100000; i++){
//         // file2 << i << "   " << 1.48e-2 << endl;          // flexion
//          file2 << i << "   " << 666e-6 << endl;             // traction
//     }
//     auto init = chrono::high_resolution_clock::now();

//     /////////////////////////////////////////////     Cas flexion 2D    /////////////////////////////////////////////

//     // double length = 2; //120; //2.0; //5; // Longueur totale de la poutre (m)
//     // double E = 30e6; //200e9;  // 72e9;  // Module de Young (en Gpa)
//     // double I = 3.6e-6; //200; //3.6e-6; // 1e-8/12;//3.6e-6;  //4e-5; // Moment d'inertie (m^4)
//     // int nbOfBeam = 1;   // Nombre d'éléments de poutre
//     // int max_iterations = 500000;
//     // double tolerance = 1e-7;


//      /////////////////////////////////////////////     Cas flexion 3D    /////////////////////////////////////////////

//     double length = 2; //120; //2.0; //5; // Longueur totale de la poutre (m)
//     double E = 200e9;  // 72e9;  // Module de Young (en Gpa)
//     double I = 3.6e-6; //200; //3.6e-6; // 1e-8/12;//3.6e-6;  //4e-5; // Moment d'inertie (m^4)
//     int nbOfBeam = 1;   // Nombre d'éléments de poutre
//     int max_iterations = 100000;
//     double tolerance = 1e-7;
//     double A = 200e-3 * 60e-3; //10;    // section en mètre carré

//     Eigen::MatrixXd coords(2,3);
//     // double Iy = 0; double Iz = 0; double G = 10000; double J = 50; 
//     double Iy = 3.6e-6; double Iz = 3.6e-6; double G = 73e9; double J = 50; 
//     coords << 0,   0,   0,
//               2,  0,   0;



//     /////////////////////////////////////////////     Cas tracion     /////////////////////////////////////////////

//     // double length = 3; //20;                 // Longueur totale de la poutre (en in)
//     // double E = 200e9;  //30e6;                    // Module de Young (en psi)
//     // double I = 0;                      // Moment d'inertie (in^4): Le moment d'inertie est nul dans le cas de la traction
//     // int nbOfBeam = 1;                   // Nombre d'éléments de poutre
//     // int max_iterations = 100000;
//     // double tolerance = 1e-7;
//     // double A = 12.5e-3 * 12.5e-3; //10;    // section en mètre carré
//     // double C = 1; 
//     // double S = 0;
//     // // double C,S;   // A revoir pour le cas où C et S dépendent des éléments poutres orientés différemment

//     // Eigen::MatrixXd coords(2,3);
//     // // double Iy = 0; double Iz = 0; double G = 10000; double J = 50; 
//     // double Iy = 100; double Iz = 100; double G = 10000; double J = 50; 
//     // // coords << 0,   0,   0,
//     // //           20,  0,   0;

//     // coords << 0,   0,   0,
//     //           3,  0,   0;
//     // // double longueur = 100;
//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//     // Eigen::MatrixXd K_global = global_stiffness_matrix_beam1D(length, E, I, nbOfBeam);                      //
//     // Eigen::MatrixXd K_global = global_stiffness_matrix_beam2D(length, A, C, E, I, S, nbOfBeam);                //
//     Eigen::MatrixXd K_global = global_stiffness_matrix_beam3D(coords, A, E, Iy, Iz, G, J, length, nbOfBeam);
//     Eigen::VectorXd  F = global_load_vector(nbOfBeam);                                                         


//     // Eigen::MatrixXd K_global_dense = global_stiffness_matrix(length, I, nbOfBeam);                   
//     // Eigen::SparseMatrix<double> K_global_sparse = K_global_dense.sparseView(); // Conversion en format sparse


//     Eigen::VectorXd U = conjugate_gradient(K_global, F, max_iterations, tolerance);                              
//     // Eigen::VectorXd U = conjugate_gradient(K_global_sparse, F, max_iterations, tolerance);

//     // Eigen::VectorXd U = gauss_seidel(K_global, F, max_iterations, tolerance);

//     // Afficher la matrice de rigidité globale
//     // std::cout << "Matrice de rigidité globale_dense :\n" << K_global_dense << std::endl;
//     // std::cout << "Matrice de rigidité globale_sparse :\n" << K_global_sparse << std::endl;
//     // std::cout << "Matrice de rigidité globale :\n" << K_global << std::endl;
//     // std::cout << "Force globale: \n" << F << std::endl;

//     // double x1 = 0; 
//     // double x2 = 2;
//     // double x2 = 120;
//     // Eigen::MatrixXd Ke = local_stiffness_matrix_beam2D(x1, x2, A, C, E, I, S);
//     // cout << "Ke: \n" << Ke/250000 << endl;

//     // Eigen::MatrixXd T(12,12);
//     // Eigen::MatrixXd Ke(12,12);
//     // Eigen::MatrixXd coords(2,3);
//     // double Iy = 100; double Iz = 100; double G = 10000; double J = 50; 
//     // coords << 0,   0,   0,
//     //           100, 0,   0;
//     // double longueur = 100;
//     // T = transformation_matrix(coords, longueur);
    
//     // Ke = local_stiffness_matrix_beam3D( coords, A, E, Iy, Iz, G, J, longueur);
//     // cout << "Transformation matrix: \n" << T << endl;
//     // cout << " Ke: \n" << Ke << endl;



//     // Afficher la solution U
//     std::cout << "Solution U :\n" << U << std::endl;                                                     

//     auto end = chrono::high_resolution_clock::now();
//     auto duration = chrono::duration_cast<chrono::seconds>(end - init).count();
//     cout << "Temps d'éxécution (en secondes): " << duration << endl;

//     return 0;
// }






