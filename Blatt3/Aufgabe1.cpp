#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;


//Funktion zum Erstellen der Hamiltonmatrix
MatrixXd Hamilton(int n, double t, double eps){
	if (n % 2 !=0)
	{
		cout << "Bitte gerade Zahl eingeben!";
	}
	MatrixXd M(n,n);
	M.setZero();
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
	return M;

}

//Funktion zur Durchführung der Jacobirotation
MatrixXd Jacobi(MatrixXd M){
	int n=M.rows();
	int z=0;
	double Offdiagonal=0;
	do
	{
		//Finde maximalen Wert der Nebendiagonale und deren Position
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

		//Erstelle Rotationsmatrix
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
		z+=1;
	} while (abs(Offdiagonal)>1e-3); //Prüfe ob Nebendiagonale gering genug

	return M;
}

//Funktion zum Durchführen des Lanczos Algorithmus
double Lanczos(MatrixXd A){

	int size = A.rows();
	vector<VectorXd> Krylov;
	VectorXd q0(size);
	q0= VectorXd::Zero(size);
	Krylov.push_back(q0);
	VectorXd k1(size);
	for (int i = 0; i < size; ++i)
	{
		k1(i)=1;
	}
	VectorXd q1(size);
	q1=k1/k1.norm();
	Krylov.push_back(q1);

	VectorXd gamma(size);
	gamma= VectorXd::Zero(size);
	gamma(0)=0;
	gamma(1)=1;

	VectorXd delta(size);
	delta= VectorXd::Zero(size);
	delta(0)=0;

	MatrixXd Eins(size,size);
	Eins.setIdentity();

	int i=1;

	do {
		delta(i)=Krylov[i].transpose()*A*Krylov[i];
		VectorXd k(size);

		k=((A-delta(i)*Eins)*Krylov[i]-gamma(i)*Krylov[i-1]);
		gamma(i+1)=k.norm();
		if (abs(gamma(i+1))<1e-7)
		{
			break;
		}
		VectorXd q(size);
		q=k*(1/gamma(i+1));
		Krylov.push_back(q);
		i+=1;
		gamma.conservativeResize(i+2);
		delta.conservativeResize(i+1);
		k.setZero();
		q.setZero();
	} while (i<size);

	//Erstelle Tridiagonalmatrix
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
	//Führe Jacobirotation mit Tridiagonalmatrix durch
	MatrixXd D;
	D = Jacobi(M);
	VectorXd x(n);
	x = D.diagonal();

	double e0=x.minCoeff();
	Krylov.clear();
	x.setZero();
	D.setZero();

	return e0;
}



int main()
{
	VectorXd e0(41);
	double epsilon = -20;
	vector<VectorXd> grund;


	for (int i = 0; i < 41; ++i)
	{
		MatrixXd Ham(50,50);
		Ham=Hamilton(50,1,epsilon);
		e0(i)=Lanczos(Ham);
		epsilon = epsilon +1;

		SelfAdjointEigenSolver<MatrixXd> es(Ham);
		VectorXd ev= es.eigenvalues();
		int minIndex;
		ev.minCoeff(&minIndex);
		VectorXd g(50);
		MatrixXd V(50,50);
		V= es.eigenvectors();
		g=V.col(minIndex);

		grund.push_back(g);
	}
	//Test ob Funktion durchläuft
	cout << e0.maxCoeff();
	cout << e0.minCoeff();

	//Speichere Resultate
	ofstream afile ("Data/one_Energie.txt", std::ofstream::out);
	afile << e0;
	afile.close();
	e0.setZero();


	ofstream bfile ("Data/one_EV.txt", std::ofstream::out);
	for (int i = 0; i < 41; ++i)
	{
		bfile << "\n" << "#epsilon:" << -20+i << "\n" << grund[i] << "\n";
	}

	//Erstelle Dichtematrix
	MatrixXd Dichte(41,50);
	for (int i = 0; i < 41; ++i)
	{
		VectorXd grundz(50);
		grundz=grund[i];
		for (int k = 0; k < 50; ++k)
		{
			Dichte(i,k)=grundz(k)*grundz(k);
		}
	}

	ofstream cfile ("Data/one_Dichte.txt", std::ofstream::out);
	cfile << Dichte;
	grund.clear();


}
