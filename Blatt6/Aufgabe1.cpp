#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <complex>


using namespace std;
using namespace Eigen;

//Funktion zur Namensgebung der txt Dateien
string give_name( const string& basename, string index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

//Eindimensionales Intervallhalbierungs Verfahren
void IH(std::function<double(VectorXd)> f, VectorXd x, VectorXd g, double &a, double &b, double &c)
{
  double d;
  while((c-a)>1e-6)
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

//Gradientenverfahren
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
        double a=-100;
        double b=0;
        double c=100;
        IH(f,x[i],gi, a,b,c);
        double lam=(a+c)/2.0;
        VectorXd xi=VectorXd::Zero(size);
        xi=x[i]+gi*lam;
        x.push_back(xi);
        i++;
        //Verhindere Endlosschleife
        if (i>100000)
        {
            cout << "Keine Konvergenz!";
            break;
        }
    } while (grad[i-1].squaredNorm()>epsilon*epsilon);
    //speichern der Ergebnisse
    ofstream afile;
    afile.open (give_name("Data/1_gradient_", name, ".txt"), ofstream::out);
    afile << "# xk(0), xk(1), f(xk), abw" << "\n";
    for (int j = 0; j < i; ++j)
    {
        afile << x[j](0) << "\t" << x[j](1) << "\t" << f(x[j]) << "\t" << (x[j]-x[i]).norm() << "\n";
    }
    afile.close();

}

//Konjugiertes Gradientenverfahren
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
        //Verhindere Endlosschleife
        if (i>10000)
        {
            cout << "Keine Konvergenz!";
            break;
        }
    } while (g_norm_squared[i]>epsilon*epsilon);
    //Speichere Ergebnisse
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
    //Implemetiere Rosenbrock Funktion und Gradient
    std::function<double(VectorXd)> rosen;
    rosen=[](VectorXd x) { return (1.0-x(0))*(1.0-x(0))+100.0*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));};
    std::function<VectorXd(VectorXd)> rosen_grad;
    rosen_grad=[](VectorXd x) {
        VectorXd grad(2);
        grad(0)=(x(0)-1.0)*2.0+400.0*x(0)*(x(0)*x(0)-x(1));
        grad(1)=200.0*(x(1)-x(0)*x(0));
        return grad;};

    //Startwert Rosenbrock
    VectorXd x0a(2);
    x0a<<-1,-1;

    //Minimierung Rosenbrock
    string name1="rosenbrock";
    Gradient(name1,rosen,rosen_grad,x0a,1e-3);
    conjugate_gradient(name1,rosen,rosen_grad,x0a,1e-3);

    //Implemetiere Funktion aus Teil b)
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
        double ex=exp(10*(x1*x2-3)*(x1*x2-3));
        double g=2*ex/((ex*(x1*x1+x2*x2)+1)*(ex*(x1*x1+x2*x2)+1));
        double g1=(10*pow(x1,3)*pow(x2,2)-30*pow(x1,2)*x2+10*x1*pow(x2,4)+x1-30*pow(x2,3));
        double g2=(10*pow(x2,3)*pow(x1,2)-30*pow(x2,2)*x1+10*x2*pow(x1,4)+x2-30*pow(x1,3));
        VectorXd grad(2);
        grad(0)=g*g1;
        grad(1)=g*g2;
        return grad;};

    //Minimierung mit Startwert 1
    VectorXd x0b1(2);
    x0b1<<1.5,2.3;
    string name2="b_1";
    Gradient(name2,f_b,f_b_grad,x0b1,1e-3);
    conjugate_gradient(name2,f_b,f_b_grad,x0b1,1e-3);

    //Minimierung mit Startwert 2
    VectorXd x0b2(2);
    x0b2<<-1.7,-1.9;
    string name3="b_2";
    Gradient(name3,f_b,f_b_grad,x0b2,1e-3);
    conjugate_gradient(name3,f_b,f_b_grad,x0b2,1e-3);

    //Minimierung mit Startwert 3
    VectorXd x0b3(2);
    x0b3<<0.5,0.6;
    string name4="b_3";
    Gradient(name4,f_b,f_b_grad,x0b3,1e-3);
    conjugate_gradient(name4,f_b,f_b_grad,x0b3,1e-3);

}
