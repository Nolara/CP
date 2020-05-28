#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <math.h>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

//Definiere Integranten für Aufgabeteil a)
double rho(double x0, double x, double y, double z){
    double rho=0;
    if (abs(x)<1 && abs(y)<1 && abs(z)<1)
    {
        rho=1/(sqrt((x0-x)*(x0-x)+y*y+z*z));
    }
    return rho;
}

//Definiere Integranten für Aufgabeteil b)
double rho2(double x0, double x, double y, double z){
    double rho=0;
    if (abs(x)<1 && abs(y)<1 && abs(z)<1)
    {
        rho=x/(sqrt((x0-x)*(x0-x)+y*y+z*z));
    }
    return rho;
}

//Funktion zur dreidimensionalen Mittelpunktsregel
double mittelpunkt(double (*f)(double,double,double,double), double a, double b, int n, double x) {
    double integral=0;

    double h=(b-a)/n;
    double x1=a;
    double x2=a+h;
    double y1=a;
    double y2=a+h;
    double z1=a;
    double z2=a+h;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                integral += h*4*f(x,(x1+x2)/2,(y1+y2)/2,(z1+z2)/2);
                x1 += h;
                x2 += h;
            }
            y1 += h;
            y2 += h;
        }
        z1 += h;
        z2 += h;
    }


    return integral;
}



int main() {

    //Aufgabeteil a) außerhalb
    VectorXd xa(70);
    VectorXd Phi(70);
    for (int i = 0; i < 70; ++i)
    {
        xa(i)=1.1+0.1*i;
        Phi(i)=mittelpunkt(rho,-1,1,1000,xa(i));
    }
    ofstream afile ("Data/two_1_out.txt", std::ofstream::out); // Erstelle txt Datei
    for (int i = 0; i < 70; ++i)
    {
        afile << xa(i) << "\t" << Phi(i) << "\n";
    }
    afile.close();

    //Aufgabeteil a) innerhalb
    VectorXd xa2(11);
    VectorXd Phi2(11);
    for (int i = 0; i < 11; ++i)
    {
        xa2(i)=0.1*i;
        Phi2(i)=mittelpunkt(rho,-1,1,1000,xa2(i));
    }
    ofstream bfile ("Data/two_1_in.txt", std::ofstream::out); // Erstelle txt Datei
    for (int i = 0; i < 11; ++i)
    {
        bfile << xa2(i) << "\t" << Phi2(i) << "\n";
    }
    bfile.close();

    //Aufgabeteil b) außerhalb
    VectorXd xb(70);
    VectorXd Phi3(70);
    for (int i = 0; i < 70; ++i)
    {
        xb(i)=1.1+0.1*i;
        Phi3(i)=mittelpunkt(rho2,-1,1,1000,xb(i));
    }
    ofstream cfile ("Data/two_2_out.txt", std::ofstream::out); // Erstelle txt Datei
    for (int i = 0; i < 70; ++i)
    {
        cfile << xb(i) << "\t" << Phi3(i) << "\n";
    }
    cfile.close();

    //Aufgabeteil b) innerhalb
    VectorXd xb2(11);
    VectorXd Phi4(11);
    for (int i = 0; i < 11; ++i)
    {
        xb2(i)=0.1*i;
        Phi4(i)=mittelpunkt(rho2,-1,1,1000,xb2(i));
    }
    ofstream dfile ("Data/two_2_in.txt", std::ofstream::out); // Erstelle txt Datei
    for (int i = 0; i < 11; ++i)
    {
        dfile << xb2(i) << "\t" << Phi4(i) << "\n";
    }
}
