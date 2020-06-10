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

void RK (string name, std::function<VectorXd(VectorXd, double m)> f, VectorXd &y0, double h, double tf, double m) {
    int n=tf/h;
	int size=y0.size();
	int dim =size/2;
    MatrixXd y=MatrixXd::Zero(size,n+1);
	y.col(0)=y0;

    VectorXd E(n);


    for (int i = 0; i < n; ++i)
    {
        E(i)=0.5*m*(y.col(i).tail(dim).squaredNorm()+y.col(i).head(dim).squaredNorm());
        VectorXd k1=h*f(y.col(i).head(dim),m)/m;
        VectorXd k2=h*f(y.col(i).head(dim)+0.5*k1,m)/m;
        VectorXd k3=h*f(y.col(i).head(dim)+0.5*k2,m)/m;
        VectorXd k4=h*f(y.col(i).head(dim)+k3,m)/m;
        VectorXd v_np1=(k1+2*k2+2*k3+k4)*1/6+y.col(i).tail(dim);
        VectorXd y_np1(6);
        y_np1.head(dim)=y.col(i).head(dim)+y.col(i).tail(dim)*h,
        y_np1.tail(dim)=v_np1;
        y.col(i+1)=y_np1;

    }
    y0=y.col(n-1);

    ofstream afile;
    afile.open (give_name("Data/1_", name, ".txt"), ofstream::out);
    afile << "# t,x1, x2, x3, v1,v2,v3" << "\n";
    for (int j = 0; j < n; ++j)
    {
        afile << h*j << "\t" << y.col(j)(0) << "\t" << y.col(j)(1) << "\t" << y.col(j)(2) << "\t" << y.col(j)(3) << "\t" << y.col(j)(4) << "\t" << y.col(j)(5) << "\t" <<  E(j) <<"\n";
    }
    afile.close();
}

int main() {
    std::function<VectorXd(VectorXd, double)> f;
    f=[](VectorXd r, double m) { return -m*r;};

    VectorXd r01(3);
    r01 <<1,0,0;
    VectorXd v01(3);
    v01 <<0,0,0;
    VectorXd y01(6);
    y01.head(3)=r01;
    y01.tail(3)=v01;


    VectorXd r02(3);
    r02 <<2,0,0;
    VectorXd v02(3);
    v02 <<3,0,0;
    VectorXd y02(6);
    y02.head(3)=r02;
    y02.tail(3)=v02;

    VectorXd r03(3);
    r03 <<2,0,0;
    VectorXd v03(3);
    v03 <<-3,0,0;
    VectorXd y03(6);
    y03.head(3)=r03;
    y03.tail(3)=v03;

    VectorXd r04(3);
    r04<<2,4,0;
    VectorXd v04(3);
    v04 <<0,-1,3;
    VectorXd y04(6);
    y04.head(3)=r04;
    y04.tail(3)=v04;

    double h=0.01;
    double tf=5*2*M_PI;
    double m=2;

    RK("1",f,y01,h,tf,m);
    RK("2",f,y02,h,tf,m);
    RK("3",f,y03,h,tf,m);
    RK("4",f,y04,h,tf,m);

    // double hb=0.1;
    // double tfb=20*2*M_PI;

    // VectorXd r0(3);
    // r0 <<1,0,0;
    // VectorXd v0(3);
    // v0 <<0,0,0;
    // VectorXd y0(6);
    // y0.head(3)=r0;
    // y0.tail(3)=v0;
	//
	//
    // ofstream bfile ("Data/1_abw.txt", std::ofstream::out); // Erstelle txt Datei
    // VectorXd abw(20);
    // int i=0;
    // do
    // {
    //     i++;
    //     hb=hb/2;
    //     VectorXd y=y0;
    //     RK("for_abw",f,y,hb,tfb,m);
    //     abw(i)=(r01-y.head(3)).norm();
    //     bfile << hb << "\t" << abw(i) << "\n";
    //     if (i>15)
    //     {
    //         break;
    //     }
    // } while (abw(i)>1e-5);

}
