#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>



#include "profiler.h"

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

MatrixXd Zeit(int max_size){
    MatrixXd Zeit(max_size, 4);
    Profiler::init( 4 );
    for (int i = 1; i < max_size; ++i)
    {
        Profiler::start( 0 );
        MatrixXd A=Random_matrix(i);
        Profiler::stop( 0 );
        Zeit(i,0)=Profiler::getTimeInNS( 0 );

        Profiler::start( 1 );
        VectorXd b=Random_vector(i);
        Profiler::stop( 1 );
        Zeit(i,1)=Profiler::getTimeInNS( 1 );

        Profiler::start( 2 );
        PartialPivLU<MatrixXd> lu(A);
        Profiler::stop( 2 );
        Zeit(i,2)=Profiler::getTimeInNS( 2 );

        Profiler::start( 3 );
        VectorXd x(i);
        x= lu.solve(b);
        Profiler::stop( 3 );
        Zeit(i,3)=Profiler::getTimeInNS( 3 );

        Profiler::resetAll();
    }
    return Zeit;
}


int main()
{
    int max_size =1000;
    MatrixXd Results(max_size,4);
    Results = Zeit(max_size);

    ofstream afile ("Data/two_data.txt", std::ofstream::out);
    afile << "# Erzeugung Zufallsmatrix, Erzeugung Zufallsvektor, LU Zerlegung, Lösung";
    afile << Results;
    afile.close();

}
