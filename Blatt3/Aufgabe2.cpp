#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <math.h>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

MatrixXd H(VectorXd m, VectorXd k){
    int n= m.size();
    MatrixXd M(n,n);

    M(0,0)=k(0)/m(0);
    M(0,1)=- k(0)/m(1);
    M(1,0)= -k(0)/m(0);
    
    for(int i=1; i<n-1; ++i){
        M(i,i)  =(k(i-1)+k(i))/m(i);
        M(i+1,i)=- k(i)/m(i);
        M(i,i+1)= -k(i)/m(i+1);
    }
    M(n-1,n-1) = k(n-2)/m(n-1);
    M(n-2,n-1) = - k(n-2)/m(n-1);
    M(n-1,n-2) =  -k(n-2)/m(n-2);


    ofstream bfile ("Data/two_matrix.txt", std::ofstream::out);
	bfile << M;
    return M;

}



int main(){

    int N = 10;
    VectorXd m(N), k(N-1) , l(N-1);
    for(int i=1; i<N; i++){
        m(i-1)=i;
        k(i-1)=N-i;
        l(i-1)=abs(5-i)+1;          // irgendwie brauchen wir l nicht???
    }
    m[N-1]=N-1;



    MatrixXd M(N, N) ;
    M=H(m,k);

    VectorXd w(N);	
	for(int i = 0; i<N; i++){
		w(i) = sqrt(abs(M.eigenvalues()[i]));
	}
	ofstream afile ("Data/two.txt", std::ofstream::out);
	afile <<  w;

    cout << "\n \n Ende Aufgabe 2 \n \n \n";

    return 0;
}
