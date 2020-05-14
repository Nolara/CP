#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;

//Funktion definieren
double Funktion1(double x){
  return exp(-x)/x;
}

double Funktion2(double x){
  if(x == 0){
    x=x+1e-7;
    return x*sin(1.0/x);
  }
  else {
    return x*sin(1.0/x);
  }


}



double Trapezregel(double (*funptr)(double), int a, int b,  double h){
  //Anzahl Teilintervalle als Integer
   int N;
   N=(b-a)/h;
   double I;
   I=0;
   double x0;
   double xi;
   x0=0.5*h*funptr(a)+0.5*h*funptr(b);
   //Interation
   int i;
   for(i=1; i<=N-1; i++){
     xi=a+i*h;
     I = I+funptr(xi);

   }
    I=I*h+x0;

   return I;

}



double Mittelpunktsregel(double (*funptr)(double),int a, int b,  double h){
  //Anzahl Teilintervalle
   int N;
   N=(b-a)/h;
   double I=0;
   double xi;
   //Iteration
   int i;
   for(i=1; i<=N-1; i++){
     xi=a-h/2+i*h;
     I=I+funptr(xi);

   }
   I=I*h;

   return I;

}



double Simpsonsregel(double (*funptr)(double), int a, int b,  double h){
  //Anzahl Teilintervalle
   int N;
   N=(b-a)/h;
   double Ii =0;
   double x0;
   double xi;
   x0=funptr(a)+funptr(b);
   //Iteration 1
   int i;
   for(i=1; i<=N-1; i++){
     xi=a+i*h;
     Ii=Ii+funptr(xi);

   }
   Ii=2*Ii;

   double Ik=0;
   double xk;
   double xk1;
   double l;
   //Iteration 2
   int k;
   for(k=1; k<=N; k++){
     xk=a+k*h;
     xk1=a+(k-1)*h;
     l=(xk+xk1)/2;
     Ik=Ik+funptr(l);

   }
   Ik=4*Ik;

   double I;
   I=(Ii+Ik+x0)*h/6;

   return I;

}


int main()
{

// Testen der Funktion 1

    //Schreibe Werte in txt Datei
    ofstream afile ("Data/3_1_1.txt", std::ofstream::out); // Erstelle txt Datei
    afile << "#Auswertung des Integrals I1 mittels Trapezregel" << "\n";
    afile << "#Integralwert,Anzahl Intervalle, Schrittweite h, rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I1;
    double h=45.5;
    I1=Trapezregel(Funktion1,1,100,h);
    double Wert1 = I1;
    double abweichung_1=0;
    double eps = 1e-4;
    //weitere Durchläufe bis relative Abweichung 10-4; h halbieren
    do{
        int N = 99/h;
        h = h/2;
        I1 = Trapezregel(Funktion1, 1, 100, h);

        abweichung_1 = 1 - I1/Wert1; // relative Abweichung
        Wert1 = I1;
        afile << I1 << "\t";
        afile << N << "\t";
        afile << h << "\t";
        afile << (abs(abweichung_1)) << "\n";

      }while(abs(abweichung_1) >= eps);

    afile.close();


    //Schreibe Werte in txt Datei
    ofstream a1file ("Data/3_1_2.txt", std::ofstream::out); // Erstelle txt Datei
    a1file << "#Auswertung des Integrals I1 mittels Mittelpunktsregel" << "\n";
    a1file << "# Integralwert,Anzahl Intervalle, Schrittweite h, rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I2;
    h=45.5;
    I2=Mittelpunktsregel(Funktion1,1,100,h);
    double Wert2 = I2;
    double abweichung_2 =0;
    //weitere Durchläufe bis relative Abweichung 10-4; h halbieren
    do{
        int N = 99/h;
        h = h/2;
        I2 = Mittelpunktsregel(Funktion1, 1, 100, h);

        abweichung_2 = 1 - I2/Wert2; // relative Abweichung
        Wert2 = I2;
        a1file << I2 << "\t";
        a1file << N << "\t";
        a1file << h << "\t";
        a1file << abs(abweichung_2) << "\n";

      }while(abs(abweichung_2) > eps);

    a1file.close();



    //Schreibe Werte in txt Datei
    ofstream a2file ("Data/3_1_3.txt", std::ofstream::out); // Erstelle txt Datei
    a2file << "#Auswertung des Integrals I1 mittels Simpsonsregel" << "\n";
    a2file << "#Integralwert, Anzahl Intervalle, Schrittweite h, Rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I3;
    h=49.5;
    I3=Simpsonsregel(Funktion1,1,100,h);
    double Wert3 = I3;
    double abweichung_3 =0;
    //weitere Durchläufe bis relative Abweichung 10-4; h halbieren
    do{
        int N = 99/h;
        h = h/2;
        I3 = Mittelpunktsregel(Funktion1, 1, 100, h);

        abweichung_3 = 1 - I3/Wert3; // relative Abweichung
        Wert3 = I3;
        a2file << I3 << "\t";
        a2file << N << "\t";
        a2file << h << "\t";
        a2file << abs(abweichung_3) << "\n";

      }while(abs(abweichung_3) > eps);

    a2file.close();


// Testen der  Funktion 2

    //Schreibe Werte in txt Datei
    ofstream bfile ("Data/3_2_1.txt", std::ofstream::out); // Erstelle txt Datei
    bfile << "#Auswertung des Integrals I2 mittels Trapezregel" << "\n";
    bfile << "#Integralwert, Schrittweite h, rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I21;
    double h1=0.5;
    I21=Trapezregel(Funktion2,0,1,h1);
    double Wert21 = I21;
    double abweichung_21=0;
    //weitere Durchläufe bis relative Abweichung 10-4; h halbieren
    do{
        h1 = h1/2;
        I21 = Trapezregel(Funktion2, 0, 1, h1);

        abweichung_21 = 1 - I21/Wert21; // relative Abweichung
        Wert21 = I21;
        bfile << I21 << "\t";
        bfile << h1 << "\t";
        bfile << (abs(abweichung_21)) << "\n";

      }while(abs(abweichung_21) >= eps);

    bfile.close();


    //Schreibe Werte in txt Datei
    ofstream b1file ("Data/3_2_2.txt", std::ofstream::out); // Erstelle txt Datei
    b1file << "#Auswertung des Integrals I2 mittels Mittelpunktsregel" << "\n";
    b1file << "# Integralwert, Schrittweite h, rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I22;
    h1=0.5;
    I22=Mittelpunktsregel(Funktion2,0,1,h1);
    double Wert22 = I22;
    double abweichung_22 =0;
    //weitere Durchläufe bis relative Abweichung 10-4; h halbieren
    do{
        h1 = h1/2;
        I22 = Mittelpunktsregel(Funktion2, 0, 1, h1);

        abweichung_22 = 1 - I22/Wert22; // relative Abweichung
        Wert22 = I22;
        b1file << I22 << "\t";
        b1file << h1 << "\t";
        b1file << abs(abweichung_22) << "\n";

      }while(abs(abweichung_22) >= eps);

    b1file.close();



    //Schreibe Werte in txt Datei
    ofstream b2file ("Data/3_2_3.txt", std::ofstream::out); // Erstelle txt Datei
    b2file << "#Auswertung des Integrals I2 mittels Simpsonsregel" << "\n";
    b2file << "#Integralwert, Schrittweite h, Rel. Abweichung" << "\n" <<  "\n";
    //Durchlauf 1
    double I23;
    h1=0.5;
    I23=Simpsonsregel(Funktion2,0,1,h1);
    double Wert23 = I23;
    double abweichung_23 =0;
    //weitere Durchläufe bis relative Abweichung 10-7; h halbieren
    do{
        h1 = h1/2;
        I23 = Mittelpunktsregel(Funktion2, 0, 1, h1);

        abweichung_23 = 1 - I23/Wert23; // relative Abweichung
        Wert23 = I23;
        b2file << I23 << "\t";
        b2file << h1 << "\t";
        b2file << abs(abweichung_23) << "\n";

      }while(abs(abweichung_23) >= eps);

    b2file.close();




      return 0;


}
