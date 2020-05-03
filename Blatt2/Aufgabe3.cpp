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

void addcontent(const std::string& name, long long& content, bool append = false) {
    std::ofstream outfile;
    if (append)
        outfile.open(name, std::ios_base::app);
    else
        outfile.open(name);
    outfile << content << "\n";
}

VectorXd Invert(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 0 );
    MatrixXd Inv(size,size);
    Inv=M.inverse();
    x= Inv*b;
    Profiler::stop( 0 );
    long long time_inv = Profiler::getTimeInNS( 0 );
    addcontent("Data/three_t_inverse.txt", time_inv ,true);
    return x;
}
VectorXd FullLu(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 1 );
    MatrixXd FullLu(size,size);
    FullPivLU<MatrixXd> flu(M);
    x= flu.solve(b);
    Profiler::stop( 1 );
    long long time_full = Profiler::getTimeInNS( 1 );
    addcontent("Data/three_t_full_LU.txt", time_full ,true);
    return x;
}
VectorXd PartLu(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 2 );
    MatrixXd PartLu(size,size);
    PartialPivLU<MatrixXd> plu(M);
    x= plu.solve(b);
    Profiler::stop( 2 );
    long long time_part = Profiler::getTimeInNS( 2 );
    addcontent("Data/three_t_partial_LU.txt", time_part ,true);
    return x;
}


int main()
{
    int max_size =1000;

    ofstream afile ("Data/three_t_inverse.txt", std::ofstream::out);
    ofstream bfile ("Data/three_t_full_LU.txt", std::ofstream::out);
    ofstream cfile ("Data/three_t_partial_LU.txt", std::ofstream::out);
    afile.close();
    bfile.close();
    cfile.close();
    Profiler::init( 3 );

    VectorXd abw1(max_size);
    VectorXd abw2(max_size);
    VectorXd abw3(max_size);

    for (int i = 1; i < max_size; ++i)
    {

        MatrixXd M(i,i);
        VectorXd b(i);
        M=Random_matrix(i);
        b=Random_vector(i);

        VectorXd x_Inv(i);
        x_Inv = Invert(i, M, b);

        VectorXd x_FullLU(i);
        x_FullLU = FullLu(i, M, b);

        VectorXd x_parLU(i);
        x_parLU = PartLu(i, M, b);

        VectorXd diff1(i);
        diff1=x_Inv-x_FullLU;
        abw1(i)=diff1.norm();

        VectorXd diff2(i);
        diff2=x_Inv-x_parLU;
        abw2(i)=diff2.norm();

        VectorXd diff3(i);
        diff3=x_parLU-x_FullLU;
        abw3(i)=diff3.norm();

        Profiler::resetAll();
    }
    ofstream dfile ("Data/three_difference.txt", std::ofstream::out);
    dfile << "# Differenzen" << "\n" ;
    for (int i = 1; i < max_size; ++i)
    {
        dfile << abw1(i) << " " << abw2(i) << " " <<  abw3(i) << "\n";
    }
    dfile.close();
}
