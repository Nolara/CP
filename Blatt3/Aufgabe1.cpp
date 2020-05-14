#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <Eigen/Eigenvalues>
#include <thread>

using namespace std;
using namespace Eigen;



MatrixXd Hamilton(int n, double t, double eps){
	if (n % 2 !=0)
	{
		cout << "Bitte gerade Zahl eingeben!";
	}
	MatrixXd M(n,n);
	M(n/2-1,n/2-1) = eps;
	for (int i = 0; i < n; ++i)
	{
		if (i < n-1)
		{
			M(i,i+1)=-t;
			M(i+1,i)=-t;
		}
		else{
			M(i,0)=-t;
			M(0,i)=-t;
		}
	}
	ofstream afile ("Data/one.txt", std::ofstream::out);
	afile << M;
	return M;
}


MatrixXd Jacobi(MatrixXd M){
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
	} while (abs(Offdiagonal)>1e-1);

	return M;
}


double Lanczos(MatrixXd A){

	int size = A.rows();

	vector<VectorXd> Krylov;
	VectorXd q0(size);
	q0= VectorXd::Zero(size);
	Krylov.push_back(q0);
	VectorXd k1(size);
	k1=VectorXd::Random(size);
	VectorXd q1(size);
	q1=k1/k1.norm();
	Krylov.push_back(q1);

	VectorXd gamma(size+4);
	gamma(0)=0;
	gamma(1)=1;

	VectorXd delta(size);
	delta(0)=0;

	int i=1;

	do {
		delta(i)=Krylov[i].transpose()*A*Krylov[i];
		VectorXd k(size);
		MatrixXd Eins(size,size);
		Eins.setIdentity();
		k=((A-delta(i)*Eins)*Krylov[i]-gamma(i)*Krylov[i-1]);
		gamma(i+1)=k.norm();
		if (abs(gamma(i+1))<1e-7)
		{
			break;
		}
		VectorXd q(size);
		q=k*(1/gamma(i+1));
		Krylov.push_back(q);
		if (Krylov[i].dot(Krylov[i+1])>1e-7)
		{
			cout << "not orthogonal!";
		}
		i+=1;
		gamma.conservativeResize(i+2);
		delta.conservativeResize(i+1);
	} while (abs(gamma(i))>1e-7);

	int m=delta.size();
	int n=m-1;
	MatrixXd M(n,n);
	M=MatrixXd::Zero(n,n);
	M(n-1,n-1)=delta(n);
	for (int i = 0; i < n-1; ++i)
	{
		M(i,i)=delta(i+1);
		M(i,i+1)=gamma(i+2);
		M(i+1,i)=gamma(i+2);
	}
	MatrixXd D;
	D = Jacobi(M);
	VectorXd x(n);
	x = D.diagonal();

	double e0=x.minCoeff();

	return e0;
}



int main()
{
	const auto processor_count = std::thread::hardware_concurrency();
	setNbThreads(processor_count);
	VectorXd eps(41);
	VectorXd e0(41);

	for (int i = 1; i < 20; ++i)
	{
		eps(i)=-20+i;
	}

	for (int i = 0; i < 1; ++i)
	{
		MatrixXd Ham(50,50);
		Ham=Hamilton(50,1,-20);
		e0(i)=Lanczos(Ham);
		cout << Ham.eigenvalues();

	}


}
