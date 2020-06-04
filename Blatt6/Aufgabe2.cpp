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


void IH(std::function<double(VectorXd)> f, VectorXd x, VectorXd g, double &a, double &b, double &c)
{
  double d;
  while((c-a)>1e-4)
  {
    if(abs(c-b)>abs(b-a))
    {
      d = (c+b)/2;
      if(f(x+g*d) < f(x+g*b))
      {
        a = b;
        b = d;
      } else {
        c = d;
      }
    } else {
      d = (b+a)/2;
      if(f(x+g*d) < f(x+g*b))
      {
        c = b;
        b = d;
      } else {
        a = d;
      }
    }
  }
}


VectorXd Gradient(std::function<double(VectorXd)> f,std::function<VectorXd(VectorXd)> g, VectorXd x0) {

    int size =x0.size();
    VectorXd gi=VectorXd::Zero(size);
    gi=-g(x0);

    double a=-10;
    double b=1;
    double c=10;
    VectorXd lambd=VectorXd::Zero(size);
    IH(f,x0,gi, a,b,c);
    double lam=(a+c)/2.0;
    VectorXd xi=VectorXd::Zero(size);
    xi=x0+gi*lam;
    return xi;

}

void BFGS(string name,std::function<double(VectorXd)> f, std::function<VectorXd(VectorXd)> g, VectorXd x0, MatrixXd C0, const double epsilon) {
    vector<MatrixXd> c;
    c.push_back(C0);

    vector<VectorXd> b;


    int size =x0.size();
    VectorXd b0=g(x0);
    VectorXd x1=Gradient(f,g,x0);
    VectorXd b1=g(x1);

    b.push_back(b0);
    b.push_back(b1);

    vector<VectorXd> x;
    x.push_back(x0);
    x.push_back(x1);
    vector<VectorXd> p;
    VectorXd p0=-C0*b0;
    p.push_back(p0);

    vector<VectorXd> s;
    vector<VectorXd> y;

    VectorXd s0=x1-x0;
    VectorXd y0=b1-b0;

    s.push_back(s0);
    y.push_back(y0);

    int k=0;
    do
    {
        double rho=1/(s[k].dot(y[k]));

        MatrixXd Eins=MatrixXd::Identity(size, size);
        MatrixXd Ckp1=(Eins-rho*s[k]*y[k].transpose())*c[k]*(Eins-rho*y[k]*s[k].transpose())+rho*s[k]*s[k].transpose();
        c.push_back(Ckp1);
        k++;
        VectorXd pk=-c[k]*b[k];
        p.push_back(pk);
        //double a=-1000;
        //double mid=1;
        //double c=1000;
        //IH(f,x[k],pk, a,mid,c);
        //double lam=(a+c)/2.0;
        VectorXd xkp1=x[k]+p[k];
        x.push_back(xkp1);
        VectorXd sk=p[k];
        s.push_back(sk);
        VectorXd bkp1=g(xkp1);
        b.push_back(bkp1);
        VectorXd yk=bkp1-b[k];
        y.push_back(yk);
        if (k>10000)
        {
            cout << "Keine Konvergenz!";
            break;
        }
    } while (b[k].squaredNorm()>epsilon*epsilon);

    ofstream afile;
    afile.open (give_name("Data/2_", name, ".txt"), ofstream::out);
    afile << "# xk(0), xk(1), f(xk), abw" << "\n";
    for (int j = 0; j < k; ++j)
    {
        afile <<j <<"\t" << x[j](0) << "\t" << x[j](1) << "\t" << f(x[j]) << "\t" << (x[j]-x[k]).norm() << "\n";
    }
    afile.close();



}


int main() {
    std::function<double(VectorXd)> f;
    f=[](VectorXd x) { return (1-x(0))*(1-x(0))+100*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0)); };
    std::function<VectorXd(VectorXd)> g;
    g=[](VectorXd x) {
        VectorXd grad(2);
        grad(0)=(x(0)-1)*2+400*x(0)*(x(0)*x(0)-x(1));
        grad(1)=200*(x(1)-x(0)*x(0));
        return grad;};
    VectorXd x0(2);
    x0 <<-1,1;


    MatrixXd H(2,2);
    H << 802,400,400,200;
    MatrixXd C1= H.inverse();
    MatrixXd C2(2,2);
    C2 << 1/802.0, 0,0, 1/200.0;
    double f0=f(x0);
    MatrixXd C3(2,2);
    C3 << f0, 0,0, f0;
    //cout << C1 << "\n" << C2 << "\n" << C3;


    const double epsilon=1e-5;
    BFGS("1",f,g,x0,C1,epsilon);
    BFGS("2",f,g,x0,C2,epsilon);
    BFGS("3",f,g,x0,C3,epsilon);

}
