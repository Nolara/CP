#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;

string give_name( const string& basename, string index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

void BFGS(string name, std::function<double(VectorXd)> f, std::function<VectorXd(VectorXd)> g, VectorXd x0,const MatrixXd CO,const double epsilon){
    VectorXd bk = g(x0);            //k index entspricht hier 0
    VectorXd xk = x0;

    VectorXd pk=-CO*bk;
    VectorXd xkk = xk+pk;     //kk Index entspriecht k+1 hier kk=1
    VectorXd bkk=g(xkk);
    VectorXd sk=xkk-xk;
    VectorXd yk = bkk-bk;

    double rhok;
    MatrixXd Ck=CO;   

    int k=0;
    ofstream file (give_name("Data/two",name,".txt"), std::ofstream::out);
    file << " #k   xmin(0)   xmin(1)  f(x)  \n\n";
    file << k << '\t' <<xk(0)<< '\t' << xk(1)   << '\t' << f(xk) <<  '\n';

    int size = x0.size();
    MatrixXd Eins= MatrixXd::Identity(size, size);

    do{
        rhok=1/(yk.dot(sk));
        
        Ck = (Eins - rhok*sk*yk.transpose())*Ck*(Eins-rhok*yk*sk.transpose())+rhok*sk*sk.transpose();
        k=k+1;

        bk=bkk;         
        xk=xkk;
    
        pk=-Ck*bk;
        xkk = xk+pk;
        bkk=g(xkk);
        sk= xkk-xk;
        yk= bkk-bk;
        if(k>10000){
            cout << "ABBRUCH; KEINE KONVERGENZ \n\n";
            break;
        }
        file << k << '\t' << xk(0) << '\t' << xk(1)  << '\t' << f(xk) << '\n';// <<'\t'<<f(xk)<< xkk-xk << '\n';
    }while(bk.squaredNorm()>epsilon);
    file.flush();
    file.close();
    cout << name << " k= " << k << "  x_min=(" << xk(0) <<  " , " << xk(0) << ")^T    f(xmin)=" << f(xk) << '\n';
}




int main()
{   
    std::function<double(VectorXd)> f;
    f=[](VectorXd x) { return (1-x(0))*(1-x(0))+100*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0)); };
    std::function<VectorXd(VectorXd)> g;
    g=[](VectorXd x) {
        VectorXd grad(2);
        grad(0)=(x(0)-1)*2+400*x(0)*(x(0)*x(0)-x(1));
        grad(1)=200*(x(1)-x(0)*x(0));
        return -grad;};
    VectorXd x0(2);
    x0(0)=-1;
    x0(1)=1;
    double epsilon = pow(10,-5);


    MatrixXd Hesse(2,2);
    Hesse(0,0)=2+400*(-x0(1)+3*x0(0)*x0(0));
    Hesse(0,1)=-400*x0(0);
    Hesse(1,0)=Hesse(0,1);
    Hesse(1,1)=200;

    MatrixXd CO_a= Hesse.inverse();
    
    MatrixXd CO_b=MatrixXd::Identity(2,2);
    CO_b(0,0)=1.0/Hesse(0,0);
    CO_b(1,1)=1.0/Hesse(1,1);

    MatrixXd CO_c= MatrixXd::Identity(2,2);
    CO_c=CO_c*f(x0);

    BFGS("_a",f, g ,x0 ,CO_a ,epsilon);
    BFGS("_b",f, g ,x0 ,CO_b ,epsilon);
    BFGS("_c",f, g ,x0 ,CO_c ,epsilon);


    //cout << Hesse << "\n\n" << CO_a << "\n\n" << CO_b << "\n\n" << CO_c << '\n';



    cout << '\n'<< "Ende Aufgabe 2" << '\n'<< '\n';
    return 0;
}




