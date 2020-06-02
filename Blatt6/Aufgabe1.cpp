#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <complex>


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


void Gradient(string name,std::function<double(VectorXd)> f,std::function<VectorXd(VectorXd)> g, VectorXd x0, double epsilon) {

    vector<VectorXd> grad;
    vector<VectorXd> x;
    x.push_back(x0);
    int size =x0.size();
    int i=0;
    do
    {
        VectorXd gi=VectorXd::Zero(size);
        gi=-g(x[i]);
        grad.push_back(gi);
        double a=-1000;
        double b=1;
        double c=1000;
        VectorXd lambd=VectorXd::Zero(size);
        //lambd= IH(f, x[i], gi, a,b,c);
        //double lam=(lambd(0)+lambd(2))/2;
        //cout << lam;
        IH(f,x[i],gi, a,b,c);
        double lam=(a+c)/2.0;
        VectorXd xi=VectorXd::Zero(size);
        xi=x[i]+gi*lam;
        x.push_back(xi);
        i++;
        if (i>100)
        {
            cout << "Keine Konvergenz!";
            break;
        }
    } while (grad[i-1].squaredNorm()>epsilon*epsilon);
    ofstream afile;
    afile.open (give_name("Data/1_gradient_", name, ".txt"), ofstream::out);
    afile << "# xk(0), xk(1), f(xk), abw" << "\n";
    for (int j = 0; j < i; ++j)
    {
        afile << x[j](0) << "\t" << x[j](1) << "\t" << f(x[j]) << "\t" << (x[j]-x[i]).norm() << "\n";
    }
    afile.close();

}

void conjugate_gradient(string name,std::function<double(VectorXd)> f,std::function<VectorXd(VectorXd)> g, VectorXd x0, double epsilon) {

    vector<VectorXd> p;
    vector<double> g_norm_squared;
    vector<VectorXd> x;
    x.push_back(x0);
    int size =x0.size();
    VectorXd g0=VectorXd::Zero(size);
    g0=-g(x0);
    double g0_norm_squared=g0.squaredNorm();
    g_norm_squared.push_back(g0_norm_squared);
    VectorXd p0=g0;
    p.push_back(p0);
    int i=0;
    do
    {
        double a=-1000;
        double b=1;
        double c=1000;
        IH(f,x[i],p[i], a,b,c);
        double lam=(a+c)/2.0;
        VectorXd xi=VectorXd::Zero(size);
        xi=x[i]+p[i]*lam;
        x.push_back(xi);
        VectorXd gi=VectorXd::Zero(size);
        gi=-g(x[i+1]);
        double gnorm =gi.squaredNorm();
        g_norm_squared.push_back(gnorm);
        double mu=g_norm_squared[i+1]/g_norm_squared[i];
        VectorXd pi=VectorXd::Zero(size);
        pi=gi+p[i]*mu;
        p.push_back(pi);
        i++;
        if (i>100)
        {
            cout << "Keine Konvergenz!";
            break;
        }
    } while (g_norm_squared[i]>epsilon*epsilon);
    ofstream afile;
    afile.open (give_name("Data/1_conjugate_", name, ".txt"), ofstream::out);
    afile << "# xk(0), xk(1), f(xk), abw" << "\n";
    for (int j = 0; j < i; ++j)
    {
        afile << x[j](0) << "\t" << x[j](1) << "\t" << f(x[j]) << "\t" << (x[j]-x[i]).norm() << "\n";
    }
    afile.close();

}



int main() {

    std::function<double(VectorXd)> rosen;
    rosen=[](VectorXd x) { return (1-x(0))*(1-x(0))+100*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));};
    std::function<VectorXd(VectorXd)> rosen_grad;
    rosen_grad=[](VectorXd x) {
        VectorXd grad(2);
        grad(0)=(x(0)-1)*2+400*x(0)*(x(0)*x(0)-x(1));
        grad(1)=200*(x(1)-x(0)*x(0));
        return grad;};

    VectorXd x0a(2);
    x0a<<-1,1;
    string name1="rosenbrock";
    Gradient(name1,rosen,rosen_grad,x0a,1e-5);
    conjugate_gradient(name1,rosen,rosen_grad,x0a,1e-5);


    std::function<double(VectorXd)> f_b;
    f_b=[](VectorXd x) { double x1=x(0);
        double x2=x(1);
        double g=1+exp(-10*(x1*x2-3)*(x1*x2-3))/(x1*x1+x2*x2);
        double func=1/g;
        return func;};
    std::function<VectorXd(VectorXd)> f_b_grad;
    f_b_grad=[](VectorXd x) {
        double x1=x(0);
        double x2=x(1);
        double g=1+exp(-10*(x1*x2-3)*(x1*x2-3))/(x1*x1+x2*x2);
        double g_inv=-exp(-10*(x1*x2-3)*(x1*x2-3))/(g*g*(x1*x1+x2*x2));
        double inx1=-20*(x1*x2-3)*x2-2*x1/((x1*x1+x2*x2)*(x1*x1+x2*x2));
        double inx2=-20*(x1*x2-3)*x1-2*x2/((x1*x1+x2*x2)*(x1*x1+x2*x2));
        VectorXd grad(2);
        grad(0)=g_inv*inx1;
        grad(1)=g_inv*inx2;
        return grad;};

    VectorXd x0b1(2);
    x0b1<<1.5,2.3;
    string name2="b_1";
    Gradient(name2,f_b,f_b_grad,x0b1,1e-3);
    conjugate_gradient(name2,f_b,f_b_grad,x0b1,1e-3);

    VectorXd x0b2(2);
    x0b2<<-1.7,-1.9;
    string name3="b_2";
    Gradient(name3,f_b,f_b_grad,x0b2,1e-3);

    VectorXd x0b3(2);
    x0b3<<0.5,0.6;
    string name4="b_3";
    Gradient(name4,f_b,f_b_grad,x0b1,1e-3);
}
