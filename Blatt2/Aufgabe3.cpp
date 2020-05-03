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

MatrixXd Inverse(int max_size, MatrixXd M, VectorXd b){
    MatrixXd Zeit(max_size, 3);
    Profiler::init( 3 );
    for (int i = 1; i < max_size; ++i)
    {
        Profiler::start( 0 );
        MatrixXd A=Random_matrix(i);
        MatrixXd b=Random_vector(i);
        Profiler::stop( 0 );
        Zeit(i,0)=Profiler::getTimeInNS( 0 );

        Profiler::start( 1 );
        A.partialPivLu();
        Profiler::stop( 1 );
        Zeit(i,1)=Profiler::getTimeInNS( 1 );

        Profiler::start( 2 );
        VectorXd x(i);
        x= A.partialPivLu().solve(b);
        Profiler::stop( 2 );
        Zeit(i,2)=Profiler::getTimeInNS( 2 );
    }
    return Zeit;
}


int main()
{
    int max_size =10000;
    MatrixXd Results(max_size,3);
    Results = Zeit(max_size);

    ofstream afile ("Data/two_data.txt", std::ofstream::out);
    afile << Results;
    afile.close();

}
