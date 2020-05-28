#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <complex>


using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
const complex<double> I(0.0,1.0);
const complex<double> pi(M_PI,0.0);


dcomp fexp(double x) {
    return exp(-x*x/2);
}
dcomp frechteck(double x) {
    dcomp f;
    if (sin(x)<0) {
        f=-1;
    } else {
        f=1;
    }
    return f;
}

//Erstellen eines Vektors mit Funktionswerten
VectorXcd make_f (dcomp (*f)(double),int m, double xa, double xb) {
    double l=xb-xa;
    int n = pow(2,m);
    double delta=l/(n);

    VectorXd x(n);
    VectorXcd vec(n);
    for (int i = 0; i < n; ++i)
    {
        x(i)=xa+delta*i;
        vec(i)=f(x(i));
    }
    return vec;
}

//Zuordnung der Zahlen durch Umordnung der binären Zahlen
int reverse(int n, int m) {
    VectorXd binaryNumber(m);
    binaryNumber.setZero();
    for (int i = 0; i < m; ++i)
    {
        binaryNumber(i) = n % 2;
        n = n / 2;
    }
    VectorXd reversed(m);
    for (int i = 0; i < m; ++i)
    {
        reversed(i)=binaryNumber(m-i-1);
    }
    int decimalNumber = 0;
    int base = 1;
    for (int i = 0; i < m; ++i)
    {
        decimalNumber += reversed(i)*base;
        base = base*2;
    }
    return decimalNumber;
}

//FFT Algorithmus
VectorXcd FFT( int m, VectorXcd f) {
    int n = pow(2,m);
    vector<MatrixXcd> s;
    MatrixXcd s0(n,n);
    for (int j = 0; j < n; ++j)
    {
        for (int l = 0; l < n; ++l)
        {
            int lquer =reverse(l,m);
            s0(j,l)=f(lquer);
        }
    }
    s.push_back(s0);
    for (int k = 1; k <= m; ++k)
    {
        MatrixXcd sk(n,n);
        sk.setZero();
        for (int j = 0; j < pow(2.0, k); ++j)
        {
            dcomp J=j;
            for (int l = 0; l < pow(2.0, m-k); ++l)
            {
                sk(j,l)=s[k-1](j,2*l)+exp(2.0*pi*I*J/pow(2.0, k))*s[k-1](j,2*l+1);
                for (int i = 0; i < pow(2.0, m-k); ++i)
                {
                    sk(j+i*pow(2,k),l)=sk(j,l);
                }
            }
        }
        //cout << sk <<"\n" << "\n";
        s.push_back(sk);
    }
    VectorXcd Fj(n);
    for (int j = 0; j < n; ++j)
    {
        Fj(j)=s[m](j,0);
    }

    return Fj;
}

VectorXcd sortiere(VectorXcd &f, int m)
{
    int n=pow(2,m);
    VectorXcd g = f;
    for (int i=0; i<n/2; i++)
    {
        f(i) = g(i+n/2);
        f(i+n/2) = g(i);
    }
    return f;
}

VectorXcd Phasenfaktor(VectorXcd f, double xa, double xb, int m) {
    dcomp L=xb-xa;
    int n = pow(2,m);
    dcomp N=n;
    dcomp delta=L/(N);
    VectorXcd g=f;
    dcomp xA=xa;
    for (int j = 0; j < n; ++j)
    {
        dcomp J=j;
        f(j)=g(j)*delta/(2.0*pi)*exp(-2.0*pi*I*J*xA/L);
    }
    return f;
}


int main() {

    //Verifizieren des Algorithmus
    ofstream afile ("Data/1_1.txt", std::ofstream::out);
    for (int m = 3; m < 5; ++m)
    {
        int n=pow(2,m);
        dcomp N=n;

        VectorXcd f1=VectorXcd::Zero(n);
        for (int i = 0; i < n; ++i)
        {
            f1(i)=sqrt(1+i);
        }
        VectorXcd Fdirekt = VectorXcd::Zero(n);
        for(double j=0; j<n; j++){
            dcomp J=j;
            for (double l=0; l<n; l++){
                dcomp L=l;
                Fdirekt(j) += exp(2.0*pi*I*L*J/N)*f1(l);
            }
        }
        VectorXcd F=FFT(m,f1);

        afile << "#Direkte Auswertung" << "\t" << "FFT bei m=" << m << "\n" << "\n";
        afile << "#Re" << "\t" << "Im" << "\t" << "Re" << "\t" << "Im" << "\n";

        for (int i = 0; i < n; ++i)
        {
            afile << real(Fdirekt(i)) <<  "\t" << imag(Fdirekt(i)) << "\t" << real(F(i)) <<  "\t"  << imag(Fdirekt(i)) << "\n";
        }
        afile << "\n" << "\n";
    }
    afile.close();
    int m=7;
    int n=pow(2,m);
    dcomp N=n;

    //Fouriertrafo Exponentialfunktion
    VectorXcd f1= make_f(fexp,m,-10,10);
    VectorXcd F1=FFT(m,f1);
    F1=sortiere(F1, m);
    F1=Phasenfaktor(F1,-10,10, m);

    //Fouriertrafo Rechteckschwingung
    VectorXcd f2= make_f(frechteck,m,-M_PI,M_PI);
    VectorXcd F2=FFT(m,f2);
    sortiere(F2, m);
    Phasenfaktor(F2,-M_PI,M_PI, m);

    //Daten speichern
    ofstream bfile ("Data/1_2.txt", std::ofstream::out);
    bfile << "#FFT[Exponentialfunktion]" << "\t" << "FFT[Rechteckschwingung] bei m=" << m << "\n" << "\n";
    bfile << "# k"  << "\t" << "Abs" << "\t" << "Re"<< "\t" << "Im" << "k"  << "\t" << "Abs" << "\t" << "\t" << "Re" << "\t"  << "Im" << "\n";
    VectorXcd Fdirekt2 = VectorXcd::Zero(n);
    for(double j=0; j<n; j++){
        dcomp J=j;
        for (double l=0; l<n; l++){
            dcomp L=l;
            Fdirekt2(j) += exp(-2.0*pi*I*L*J/N)*f1(l);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        double kexp=(i-n/2.0)*2*M_PI/(20);
        bfile << kexp  << "\t" << abs(F1(i)) << "\t" << real(F1(i)) << "\t" << imag(F1(i)) << "\t" << abs(Fdirekt2(i)) << "\t" << real(F2(i)) << "\t" << imag(F2(i)) << "\n";
    }


}
