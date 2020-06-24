#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include "profiler.h"

using namespace std;
using namespace Eigen;

//Funktion zur Kraft
//c wird zum Zählen der Funktionsaufrufe verwendet
VectorXd f(VectorXd x, VectorXd v, double a, int &c) {
	c +=1;
	return -x-a*v;

}

//Funktion zur Namensvergebung in txt Dateien
string give_name( const string& basename, string index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

//Runge-Kutta-Verfahren für n Schritte, Svhrittweiten h
MatrixXd RK (VectorXd (*f)(VectorXd, VectorXd, double, int &), VectorXd y0, double h, double a, int  n, int &count) {
	int size=y0.size();
    int dim=size/2; //Dimension bestimmen
    MatrixXd y=MatrixXd::Zero(size,n+1); //Matrix zum Speichern von x und v
    y.col(0) =y0;

    for (int i = 0; i < n; ++i)
    {
        VectorXd k1=h*f(y.col(i).head(dim),y.col(i).tail(dim),a,count);
        VectorXd k2=h*f(y.col(i).head(dim)+0.5*k1,y.col(i).tail(dim),a,count);
        VectorXd k3=h*f(y.col(i).head(dim)+0.5*k2,y.col(i).tail(dim),a,count);
        VectorXd k4=h*f(y.col(i).head(dim)+k3,y.col(i).tail(dim),a,count);
        VectorXd v_np1=(k1+2*k2+2*k3+k4)*1/6+y.col(i).tail(dim);
        VectorXd y_np1(size);
		//Neuer Vektor y aus x(Head) und v(Tail)
        y_np1.head(dim)=y.col(i).head(dim)+y.col(i).tail(dim)*h;
        y_np1.tail(dim)=v_np1;
		//Übergeben vom Vektor y an die Matrix
        y.col(i+1)=y_np1;

    }
	//Ausgeben der Matrix y mit den ersten n Schritten
	return y;

}

//Adams-Bashforth Verfahren
void Adams_Bashforth (string name,VectorXd (*f)(VectorXd, VectorXd, double, int &), VectorXd y0, double h, double tf, double a, int &count) {
    int size=y0.size();
    int dim=size/2; //Dimension bestimmen
    int n=tf/h; //Anzahl der Schritte bestimmen
    MatrixXd y=MatrixXd::Zero(size,n+4); //Matrix für die Ergebnisse
    MatrixXd F=MatrixXd::Zero(size,n+4); //Matrix für die Funktionswerte
    MatrixXd y_start=RK(f,y0,h,a,3,count); //Bestimmung der ersten 3 Schritte mit Runge-Kutta
	//Die ersten Schritte übergeben
    for (int i = 0; i < 4; ++i)
    {
        y.col(i)=y_start.col(i);
        F.col(i).head(dim)=y.col(i).tail(dim);
        F.col(i).tail(dim)=f(y.col(i).head(dim),y.col(i).tail(dim),a,count);
    }
	//Eigentliches Runge-Kutta Verfahren für Schritte ab n=4
    for (int i = 3; i < n; ++i)
    {
        F.col(i).head(dim)=y.col(i).tail(dim);
        F.col(i).tail(dim)=f(y.col(i).head(dim),y.col(i).tail(dim),a,count);
        y.col(i+1)=y.col(i)+(55*F.col(i)-59*F.col(i-1)+37*F.col(i-2)-9*F.col(i-3))*h/24.0;
    }

	//Abspeichern der Ergebnisse

	//Zur Bestimmung der Laufzeit den folgenden Abschnitt bitte auskommentieren!
    ofstream afile;
    afile.open (give_name("Data/2_", name, ".txt"), ofstream::out);
    afile << "# t,x,v, E_kin, E_pot" << "\n";
    for (int j = 0; j < n; ++j)
    {
		double Ekin =0.5*y.col(j).tail(dim).dot(y.col(j).tail(dim));
		double Epot=a* y.col(j).head(dim).dot(y.col(j).tail(dim))+0.5*y.col(j).head(dim).dot(y.col(j).head(dim));
        afile << h*j << "\t" << y.col(j).head(dim) << "\t" << y.col(j).tail(dim) << "\t" << Ekin << "\t" << Epot << "\n";
    }
    afile.close();
}

int main()
{
	//Startvektor bestimmen
    VectorXd r01(1);
    r01 <<1;
    VectorXd v01(1);
    v01 <<0;
    VectorXd y01(2);
    y01.head(1)=r01;
    y01.tail(1)=v01;

	//Aufgabenteil a

	//Schrittweite und Zeit festlegen
    double h=0.0001;
    double tf=5*2*M_PI;

	//Werte für alpha festlegen
    double a1=0.5;
	double a2=2.0;
	double a3=4.0;
	double a4=0.0;
	double a5=-0.1;

	int count =0;

    Adams_Bashforth("1",f,y01,h,tf,a1,count);
	Adams_Bashforth("2",f,y01,h,tf,a2,count);
	Adams_Bashforth("3",f,y01,h,tf,a3,count);
	Adams_Bashforth("4",f,y01,h,tf,a4,count);
	Adams_Bashforth("5",f,y01,h,tf,a5,count);

	//Aufgabenteil b

	//Schrittweite, Zeit  und alpha festlegen
	double h_b=0.0001;
	double t_b=20;
	double a_b=0.1;

	Adams_Bashforth("b",f,y01,h_b,t_b,a_b,count);

	//Aufgabenteil c

   //Schrittweite und alpha festlegen

	double h_c=0.0001;
	double a_c=0.1;

	//Zeiten festlegen
	double t_c_1=10;
	double t_c_2=50;
	double t_c_3=100;
	double t_c_4=1000;

	//Startwerte zum Zählen der Funktionsaufrufe
	int count_1_ab=0;
	int count_1_rk=0;
	int count_2_ab=0;
	int count_2_rk=0;
	int count_3_ab=0;
	int count_3_rk=0;
	int count_4_ab=0;
	int count_4_rk=0;

	//Profiler zur Zeitmessung initialisieren
	Profiler::init( 8 );

	Profiler::start( 0 );
	Adams_Bashforth("c_1",f,y01,h_c,t_c_1,a_c,count_1_ab);
	Profiler::stop( 0 );

	Profiler::start( 1 );
	RK(f,y01,h_c,a_c, t_c_1/h_c,count_1_rk);
	Profiler::stop( 1 );

	Profiler::start( 2 );
	Adams_Bashforth("c_2",f,y01,h_c,t_c_2,a_c,count_2_ab);
	Profiler::stop( 2 );

	Profiler::start( 3 );
	RK(f,y01,h_c,a_c, t_c_2/h_c,count_2_rk);
	Profiler::stop( 3 );

	Profiler::start( 4 );
	Adams_Bashforth("c_3",f,y01,h_c,t_c_3,a_c,count_3_ab);
	Profiler::stop( 4 );

	Profiler::start( 5 );
	RK(f,y01,h_c,a_c, t_c_3/h_c,count_3_rk);
	Profiler::stop( 5 );

	Profiler::start( 6 );
	Adams_Bashforth("c_4",f,y01,h_c,t_c_4,a_c,count_4_ab);
	Profiler::stop( 6 );

	Profiler::start( 7 );
	RK(f,y01,h_c,a_c, t_c_4/h_c,count_4_rk);
	Profiler::stop( 7 );

	//Werte der Laufzeit und Funktionsaufrufe abspeichern
	ofstream dfile ("Data/2_Vergleich.txt", std::ofstream::out);
    dfile << "# Zeit, Laufzeit AB, Funktionsaufrufe AB, Laufzeit RK, Funktionsaufrufe RK, " << "\n" ;
    dfile << t_c_1 << "\t" << Profiler::getTimeInS( 0 ) << "\t" << count_1_ab << "\t" << Profiler::getTimeInS( 1 ) << "\t" << count_1_rk << "\n";
	dfile << t_c_2 << "\t" << Profiler::getTimeInS( 2 ) << "\t" << count_2_ab << "\t" << Profiler::getTimeInS( 3 ) << "\t" << count_2_rk << "\n";
	dfile << t_c_3 << "\t" << Profiler::getTimeInS( 4 ) << "\t" << count_3_ab << "\t" << Profiler::getTimeInS( 5 ) << "\t" << count_3_rk << "\n";
	dfile << t_c_4 << "\t" << Profiler::getTimeInS( 6 ) << "\t" << count_4_ab << "\t" << Profiler::getTimeInS( 7 ) << "\t" << count_4_rk << "\n";
    dfile.close();


}
