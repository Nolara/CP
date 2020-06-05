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

MatrixXd RK (std::function<VectorXd(VectorXd, VectorXd, double)> f, VectorXd y0, double h, double a) {
    int size=y0.size();
    int dim=size/2;
    MatrixXd y=MatrixXd::Zero(size,4);
    y.col(0) =y0;

    for (int i = 0; i < 3; ++i)
    {
        VectorXd k1=h*f(y.col(i).head(dim),y.col(i).tail(dim),a);
        VectorXd k2=h*f(y.col(i).head(dim)+0.5*k1,y.col(i).tail(dim),a);
        VectorXd k3=h*f(y.col(i).head(dim)+0.5*k2,y.col(i).tail(dim),a);
        VectorXd k4=h*f(y.col(i).head(dim)+k3,y.col(i).tail(dim),a);
        VectorXd v_np1=(k1+2*k2+2*k3+k4)*1/6+y.col(i).tail(dim);
        VectorXd y_np1(size);
        y_np1.head(dim)=y.col(i).head(dim)+y.col(i).tail(dim)*h,
        y_np1.tail(dim)=v_np1;
        y.col(i+1)=y_np1;

    }
    return y;
}

void Adams_Bashforth (string name,std::function<VectorXd(VectorXd, VectorXd, double)> f, VectorXd y0, double h, double tf, double a) {
    int size=y0.size();
    int dim=size/2;
    int n=tf/h;
    MatrixXd y=MatrixXd::Zero(size,n+4);
    MatrixXd F=MatrixXd::Zero(size,n+4);
    MatrixXd y_start=RK(f,y0,h,a);
    for (int i = 0; i < 4; ++i)
    {
        y.col(i)=y_start.col(i);
        F.col(i).head(dim)=y.col(i).tail(dim);
        F.col(i).tail(dim)=f(y.col(i).head(dim),y.col(i).tail(dim),a);
    }
    for (int i = 3; i < n; ++i)
    {
        F.col(i).head(dim)=y.col(i).tail(dim);
        F.col(i).tail(dim)=f(y.col(i).head(dim),y.col(i).tail(dim),a);
        y.col(i+1)=y.col(i)+(55*F.col(i)-59*F.col(i-1)+37*F.col(i-2)-9*F.col(i-3))*h/24.0;
    }
    ofstream afile;
    afile.open (give_name("Data/2_", name, ".txt"), ofstream::out);
    afile << "# t,x,v" << "\n";
    for (int j = 0; j < n; ++j)
    {
        afile << h*j << "\t" << y.col(j)(0) << "\t" << y.col(j)(1) << "\n";
    }
    afile.close();
}

int main()
{
    std::function<VectorXd(VectorXd, VectorXd, double)> f;
    f=[](VectorXd x,VectorXd v, double a) { return -x-a*v;};

    VectorXd r01(1);
    r01 <<1;
    VectorXd v01(1);
    v01 <<0;
    VectorXd y01(2);
    y01.head(1)=r01;
    y01.tail(1)=v01;
    double h=0.0001;
    double tf=5*2*M_PI;
    double a=0.1;

    Adams_Bashforth("1",f,y01,h,tf,a);
}
