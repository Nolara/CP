#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

double trapez(double (*f)(double), double a, double b, int n) {
    double integral;
    integral=0;
    double h;
    h=(b-a)/n;
    double x1;
    double x2;
    x1=a;
    x2=a+h;
    if (isnan(f(x1)))
    {
        x1+=1e-10;
    }
    for (int i = 0; i < n; ++i)
    {
        integral += f(x1)+f(x2);
        x1 += h;
        x2 += h;
    }
    integral=integral*h/2.0 ;
    return integral;
}

//Definiere Funktion zur Mittelpunktsregel
double mittelpunkt(double (*f)(double), double a, double b, int n) {
    double integral;
    integral=0;
    double h;
    h=(b-a)/n;
    double x1;
    double x2;
    x1=a;
    x2=a+h;
    if (isnan(f(x1)))
    {
        x1+=1e-10;
    }
    for (int i = 0; i < n; ++i)
    {
        integral += f((x1+x2)/2.0);
        x1 += h;
        x2 += h;
    }
    integral=integral*h;
    return integral;
}

//Definiere Funktion zur Simpsonregel
double simpson(double (*f)(double), double a, double b, int n) {
      double integral;
      integral=0;
      double h;
      h=(b-a)/n;
      double x1;
      double x2;
      x1=a;
      x2=a+h;
      if (isnan(f(x1)))
      {
          x1+=1e-10;
      }
      for (int i = 0; i < n; ++i)
      {
          integral += f(x1)+f(x2)+4*f((x1+x2)/2.0);
          x1 += h;
          x2 += h;
      }
      integral=integral*h/6.0 ;
      return integral;
}

double Funktion1(double x) {
    return exp(x)/x;
}
double PartI2(double x) {
    return 1+x*x/6;
}
double PartI3(double x) {
    return (exp(x)-exp(-x))/x;
}
double PartII(double x) {
    return exp(-x)/sqrt(x);
}
double PartII2(double x) {
    return exp(-1/x)/(sqrt(x)*x);
}
double PartIII(double x) {
    return sin(x)/x;
}
double PartIII2(double x) {
    return (sin(x)-sin(-x))/x;
}

double calculate(double (*func)(double), double (*integrater)(double (*f)(double), double a, double b, int n), double unten, double oben, int n, double ex){
		VectorXd integral(2);
		int i=0;
		integral(0)=0.0;
		do
		{
			i+=1;
			integral(i)=integrater(func , unten ,oben , n);
			n=2*n;
			integral.conservativeResize(i+2);

		} while (abs((integral(i)-integral(i-1))/integral(i-1))>ex);
		return integral(i);

    }


int main(){

	const int n=1e3;


  double delta=1e-6;
  double delta_z=0-delta;
  double integral1 = calculate(Funktion1,simpson,-1,-delta,n,1e-7)+calculate(Funktion1,simpson,delta,1,n,1e-7)+calculate(Funktion1,mittelpunkt,delta_z,0,n,1e-7);


	const double xmax=1e3;
	double integral21 = calculate(PartII,mittelpunkt,0, 1, n, 1e-5)+calculate(PartII2,mittelpunkt,0, 1, n, 1e-7) ;
  double integral22 = calculate(PartII,mittelpunkt,0, 1, n, 1e-5)+calculate(PartII,mittelpunkt,1, xmax, n, 1e-5);


	double integral3=2*calculate(PartIII2,mittelpunkt,0, 1, n, 1e-7)+calculate(PartIII,mittelpunkt,1, xmax, n, 1e-7);



	ofstream afile ("Data/one.txt", std::ofstream::out); // Erstelle txt Datei
  afile << integral1 << "\n";
	afile << integral21 << "\n";
	afile << integral22 << "\n";
	afile << integral3 << "\n";

}
