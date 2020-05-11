#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <string>
#include <sstream>

using namespace std;
using namespace Eigen;


//Definiere Funktionen zum Integrieren
double PartI(double x) {
    return exp(-x)/x;
}

double PartII(double x) {
    return x * sin(1/x);
}

//Definiere Funktion zur Trapezregel
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

void calculate(double (*func)(double), double (*integrater)(double (*f)(double), double a, double b, int n), double unten, double oben, int n){
    double t1=integrater(func, unten, oben, n);
    double t2=integrater(func, unten, oben, 2*n);
    while (abs((t1-t2)/t1)>0.0001)
    {
        n=n*2;
        t1=integrater(func, unten, oben, n);
        t2=integrater(func, unten, oben, 2*n);
    }
    string name;
    if (integrater==simpson) {
        name="Simpsonregel";
    } else if (integrater==trapez) {
        name="Trapezregel";
    } else {
        name="Mittelpunktsregel";
    }

    string namepart;
    if (func==PartI)
    {
        namepart= "i";
    }
    else
    {
        namepart= "ii";
    }

    double h=(oben-unten)/n;
    
    std::ofstream afile;
    afile.open("Data/three.txt", std::ios_base::app);
    afile << "Ergebniss der " << name << " aus Aufgabenteil " << namepart << " ist:" << "\n";
    afile << t1 <<"\n";
    afile << "Dazu wurde das Intervall in " << n << " Stücke der Breite " << h << " unterteilt."<< "\n" << "\n";
}



int main() {

    double i_unten=1;
    double i_oben=100;

    double ii_unten=0;
    double ii_oben=1;

    int n_start=2;


    ofstream afile ("Data/three.txt", std::ofstream::out);

    calculate(PartI,trapez, i_unten, i_oben, n_start);
    calculate(PartI,mittelpunkt, i_unten, i_oben, n_start);
    calculate(PartI,simpson, i_unten, i_oben, n_start);
    calculate(PartII,trapez, ii_unten, ii_oben, n_start);
    calculate(PartII,mittelpunkt, ii_unten, ii_oben, n_start);
    calculate(PartII,simpson, ii_unten, ii_oben, n_start);

}
