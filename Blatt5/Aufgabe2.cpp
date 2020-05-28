#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

double f(double x) {
    return x*x-2;
}

double Ableitung1 (double (*f)(double), double x) {
    double delta=1e-9;
    double abl = (f(x+delta)-f(x-delta))/(2*delta);
    return abl;
}

double Ableitung2 (double (*f)(double), double x) {
    double delta=1e-3;
    double abl=(f(x+delta)+f(x-delta)-2*f(x))/(delta*delta);
    return abl;
}
void IH (double (*f)(double), double a, double b, double c) {
    int n=0;
    int m=0;
    do
    {
        n+=1;
        if (b-a > c-b) {
            double d=(a+b)/2.0;
            if (f(d)<f(b))
            {
                m+=2;
                a=d;
            }
            else
            {
                c=b;
                b=d;

            }
        } else {
            double d=(b+c)/2.0;
            if (f(d)<f(b)) {
                m+=2;
                a=b;
                b=d;
            } else {
                c=d;
            }
        }
    } while (c-a>1e-9);
    ofstream bfile ("Data/2_Intervallhalbierung.txt", std::ofstream::out); // Erstelle txt Datei
    bfile << "Anzahl Iterationsschritte:" << n << "\n";
    bfile << "Minimum bei x=" << (c-a)/2 << " mit Funktionswert f(x)=" << f((c-a)/2) << "\n";
    bfile << "Anzahl Funktionsaufrufe:" << m << "\n";
    bfile.close();
}


void Newton (double (*f)(double), double x0) {
    int n =0;
    int m=0;
    double delta_x=0;

    delta_x=-Ableitung1(f,x0)/Ableitung2(f,x0);
    do
    {
        n +=1;

        delta_x=-Ableitung1(f,x0)/Ableitung2(f,x0);
        m+=5;
        if (f(x0+delta_x)>f(x0))
        {
            cout << "Keine Konvergenz!" << "\n";
            break;
        }
        m+=2;
        x0+=delta_x;
    } while (abs(delta_x)>1e-9);
    ofstream afile ("Data/2_Newton.txt", std::ofstream::out); // Erstelle txt Datei
    afile << "Anzahl Iterationsschritte:" << n << "\n";
    afile << "Minimum bei x=" << x0 << " mit Funktionswert f(x)=" << f(x0) << "\n";
    afile << "Anzahl Funktionsaufrufe:" << m << "\n";
    afile.close();
}

int main() {

    Newton(f, 1);
    IH(f, -0.5,-0.1,2);
}
