#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>


using namespace std;
using namespace Eigen;

string give_name( const string& basename, string index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}


VectorXd funkt(VectorXd r, double m){
    return -m*r;
}

void RungeKutta(string name,std::function<VectorXd(VectorXd,double m)> f,VectorXd &r0,VectorXd &v0, double &h,double t_stopp,double m){
    ofstream file (give_name("Data/one",name,".txt"), std::ofstream::out);
    file << "#n  y, y2 y3 y v1 v2 v3 v Eges \n\n";

    VectorXd k1,k2,k3,k4,yn,vn;
    double tn;
    double E;
    yn=r0;
    vn=v0;
    E=0.5*m*(vn.squaredNorm()+yn.squaredNorm());
    int N=t_stopp/h;
    file << 0 << '\t' << yn.transpose() <<'\t' << yn.norm() <<'\t' << vn.transpose() <<'\t' << vn.norm()<<'\t'<<E<<'\n';
    for(int n=1; n<=N; n++){
        tn = n*h;
        k1=h*f(yn,m)/m;
        k2=h*f(yn+k1/2,m)/m;
        k3=h*f(yn+k2/2,m)/m;
        k4=h*f(yn+k3,m)/m;
        yn= yn + vn*h;
        vn= vn +(k1+2*k2+2*k3+k4)/6;
        E=0.5*m*(vn.squaredNorm()+yn.squaredNorm());
        file << tn << '\t' << yn.transpose() <<'\t' << yn.norm()<<'\t' << vn.transpose() <<'\t' << vn.norm() <<'\t'<<E << '\n';
    }
    file.flush();
    file.close();
}


// Hier wie in void RungeKutta, nur ohne Abspeicherung der Daten und ohne tn,E -> Schnelleres rechnenfür die b)
void RungeKutta2(std::function<VectorXd(VectorXd,double m)> f,VectorXd &r0,VectorXd &v0, double &h,double t_stopp,double m){
    VectorXd k1,k2,k3,k4,yn,vn;
    yn=r0;
    vn=v0;
    int N=t_stopp/h;
    for(int n=1; n<=N; n++){
        k1=h*f(yn,m)/m;
        k2=h*f(yn+k1/2,m)/m;
        k3=h*f(yn+k2/2,m)/m;
        k4=h*f(yn+k3,m)/m;
        yn= yn + vn*h;
        vn= vn +(k1+2*k2+2*k3+k4)/6;
    }
    r0=yn;  // Um die Abweichung r0-ri außerhalb dieser Fkt bestimmen zu können
}


int main()
{
    double h = 1./100;
    double t_stopp = 20.;
    double m = 1.;

    VectorXd v0(3);
    v0<< 0,0,0;
    VectorXd r0 = VectorXd::Random(3);
    RungeKutta("1",funkt, r0,v0, h,t_stopp,m);      //r0 beliebig, v0=0
    double h1 = 1./10;
    RungeKutta("1b",funkt, r0,v0, h1,t_stopp,m);
    double h2 = 1./1000;
    RungeKutta("1c",funkt, r0,v0, h2,t_stopp,m);


    v0 = VectorXd::Random(3);
    r0 = v0;
    RungeKutta("2",funkt, r0,v0, h,t_stopp,m);  //v0 beliebig, r0=v0

    v0 = VectorXd::Random(3);
    r0 = VectorXd::Random(3);
    RungeKutta("3",funkt, r0,v0, h,t_stopp,m);   // r0, v0 beliebig


    /// b)
    ofstream fileb ("Data/one_abw.txt", std::ofstream::out);
    fileb << "# h  abw" << "\n\n";
    t_stopp = 20*M_PI;
    v0<< 0,0,0;
    r0 = VectorXd::Random(3);
    VectorXd ri = r0;
    double abw;
    h = 1./100;
    do{
        RungeKutta2(funkt, ri,v0,h,t_stopp,m);
        abw=(r0-ri).norm();
        fileb << h << '\t' << abw << '\n';
        ri=r0;
        h=h/2;
    }while(abw>pow(10,-2)); //Hier muss eigentlich pow(10,-5) aufgrund langer Laufzeit zeitweise geändert
    fileb.flush();
    fileb.close();



    cout << '\n'<< "Ende Aufgabe 1" << '\n'<< '\n';
}
