#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

using namespace std;
using namespace Eigen;

MatrixXd Random_matrix(int n){
    MatrixXd A(n,n);
    A= MatrixXd::Random(n,n);
    return A; // Ausgabe des Ergebnis-Vektors
}

VectorXd Random_vector(int n){
    VectorXd v(n);
    v=VectorXd::Random(n);
    return v; // Ausgabe des Ergebnis-Vektors
}


int main()
{
    int n=5;
    MatrixXd A1=Random_matrix(n);
    MatrixXd b=Random_vector(n);


    A1.partialPivLu();
    VectorXd x(n);
    x= A1.partialPivLu().solve(b);

    ofstream afile ("Data/two_data.txt", std::ofstream::out);

    afile << "# Matrix A:" << "\n";

    for (int i = 0; i < n; ++i) //Schreibe Matrix ATA in txt Datei
    {
        for (int j = 0; j < n; ++j)
        {
            afile << A1(i,j) << "\t";
        }
        afile << "\n";
    }

    afile  << "\n" << "# Vector b:" << "\n";

    for(int i=0; i<n; i++){
        afile << b(i) << "\n"; //Schreibe ATy in txt Datei
    }

    afile  << "\n" << "# Vector x:" << "\n";

    for(int i=0; i<n; i++){
        afile << x(i) << "\n"; //Schreibe ATy in txt Datei
    }

    afile.close();

}
