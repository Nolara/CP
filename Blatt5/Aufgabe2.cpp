#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

double f(double x, int &m) {
    m +=1;
    return x*x-2;
}

//Erste Ableitung
double Ableitung1 (double (*f)(double, int &), double x, int &m) {
    double delta=1e-9;
    double abl = (f(x+delta,m)-f(x-delta,m))/(2*delta);
    return abl;
}

//Zweite Ableitung
double Ableitung2 (double (*f)(double,int &), double x, int &m) {
    double delta=1e-3;
    double abl=(f(x+delta,m)+f(x-delta,m)-2*f(x,m))/(delta*delta);
    return abl;
}

//Intervallhalbierungs-Verfahren
void IH (double (*f)(double,int &), double a, double b, double c, double accuracy) {
    int n=0;
    int m=0;
    do
    {
        n+=1;
        if (b-a > c-b) {
            double d=(a+b)/2.0;
            if (f(d,m)<f(b,m))
            {
                a=d;
            }
            else
            {
                c=b;
                b=d;

            }
        } else {
            double d=(b+c)/2.0;
            if (f(d,m)<f(b,m)) {
                a=b;
                b=d;
            } else {
                c=d;
            }
        }
    } while (c-a>accuracy);
    ofstream bfile ("Data/2_Intervallhalbierung.txt", std::ofstream::out); // Erstelle txt Datei
    bfile << "Anzahl Iterationsschritte:" << n << "\n";
    bfile << "Minimum bei x=" << (c-a)/2 << " mit Funktionswert f(x)=" << f((c-a)/2,m) << "\n";
    bfile << "Anzahl Funktionsaufrufe:" << m << "\n";
    bfile.close();
}

//Newton-Verfahren
void Newton (double (*f)(double, int &), double x0, double accuracy) {
    int n =0;
    int m=0;
    double delta_x=0;
    do
    {
        n +=1;

        delta_x=-Ableitung1(f,x0,m)/Ableitung2(f,x0,m);
        if (f(x0+delta_x,m)>f(x0,m))
        {
            cout << "Keine Konvergenz!" << "\n";
            break;
        }
        x0+=delta_x;
    } while (abs(delta_x)>accuracy);
    ofstream afile ("Data/2_Newton.txt", std::ofstream::out); // Erstelle txt Datei
    afile << "Anzahl Iterationsschritte:" << n << "\n";
    afile << "Minimum bei x=" << x0 << " mit Funktionswert f(x)=" << f(x0,m) << "\n";
    afile << "Anzahl Funktionsaufrufe:" << m << "\n";
    afile.close();
}

int main() {
    double accuracy=1e-9;

    Newton(f, 1, accuracy);
    IH(f, -0.5,-0.1,2, accuracy);
}
