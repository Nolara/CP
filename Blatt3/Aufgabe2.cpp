#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <math.h>

using namespace std;
using namespace Eigen;

MatrixXd Hamilton(VectorXd m, VectorXd k){
    int n= m.size();
    MatrixXd M(n,n);

    M(0,0)=-k(0)/m(0);
    M(0,1)=k(0)/m(0);
    M(1,0)=k(0)/m(0);

    M(n-1,n-1)=-k(n-2)/m(n-1);
    M(n-2,n-1)=k(n-2)/m(n-2);
    M(n-2,n-1)=k(n-2)/m(n-1);
    for (int i = 1; i < n-1; ++i)
    {
        M(i,i)=-(k(i-1)+k(i))/m(i);
        M(i,i+1)=k(i)/m(i);
        M(i+1,i)=k(i)/m(i+1);
    }
    ofstream bfile ("Data/two_matrix.txt", std::ofstream::out);
	bfile << M;
    return M;

}

VectorXd Jacobi(MatrixXd M){
	int n=M.rows();
	vector<VectorXd> Eigenvalues;
	VectorXd ev0(n);
	ev0=VectorXd::Zero(n);
	Eigenvalues.push_back(ev0);
	int z=0;
	double Offdiagonal=0;
	do
	{
		double max = 0;
		int row = 0;
		int col = 0;
		for (int k = 0; k < n; ++k)
		{
			for (int l = 0; l < n; ++l)
			{
				if ( abs(M(k,l)) > max && l != k)
				{
					max= abs(M(k,l));
					row = k;
    				col = l;
				}
				if (l!=k)
				{
					Offdiagonal+=abs(M(k,l));
				}
			}
		}
		if (col<row)
		{
			int a=row;
			row=col;
			col=a;
		}
		double omega=(M(col,col)-M(row,row))/(2*M(row,col));
		//cout << M <<"\n" << "\n"<< row << "\n" << col << "\n" << "\n"
		double t;
		if (omega<0) {
		    t=1/(omega-sqrt(1+omega*omega));
		} else {
		    t=1/(omega+sqrt(1+omega*omega));
		}
		double c=1/(sqrt(1+t*t));
		if (c==1)
		{
			break;
		}
		double s=c*t;
		MatrixXd Rot(n,n);
		Rot.setIdentity();
		Rot(col,col)=c;
		Rot(row,row)=c;
		Rot(row,col)=s;
		Rot(col,row)=-s;
		MatrixXd Rott(n,n);
		Rott=Rot.transpose();
		//cout << Rott << "\n" << "\n" << M << "\n" << "\n" << Rot << "\n" << "\n";
		MatrixXd Mat(n,n);
		Mat.setZero();
		M=Rott*M*Rot;
		VectorXd ev(n);
		for (int j = 0; j < n; ++j)
		{
			ev(j)=M(j,j);
		}
		Eigenvalues.push_back(ev);
		z+=1;
	} while (abs(Offdiagonal)>1e-6);

    VectorXd x(n);
	x = M.diagonal();
	ofstream afile ("Data/two.txt", std::ofstream::out);
	afile << x;


	return x;
}


int main()
{
    int size = 10;

    VectorXd m(size);
    VectorXd k(size-1);

    for (int i = 0; i < size-1; ++i)
    {
        m(i)=i+1;
        k(i) = size-(i+1);
    }
    m(size-1)=size;

    MatrixXd M(size, size);
    M=Hamilton(m,k);
    VectorXd L(size);
    L=Jacobi(M);
    double x=L.minCoeff();
    double ev=sqrt(-x);

    //cout << "\n" << M.eigenvalues() << "\n"<< "\n"<< L  ;
}
