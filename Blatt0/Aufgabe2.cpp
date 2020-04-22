#include <iostream>
#include <math.h> // exp-Funktion
#include <fstream> // Ausgabe als txt
using namespace std;


void funktion_a(double xa_0, double xa_f, int n, double* ya){ // Funktion aus Aufgabenteil a
    double dxa = (xa_f-xa_0)/n;
    for(int i=0; i<=n; i++){
        ya[i] = 1/sqrt(xa_0)-1/sqrt(xa_0+1);
        xa_0=xa_0+dxa;
  }
}

void funktion_b(double xb_0, double xb_f, int n, double* yb){ // Funktion aus Aufgabenteil a
    double dxb = (xb_f-xb_0)/n;
    for(int i=0; i<=n; i++){
        yb[i] = (1-cos(xb_0))/sin(xb_0);
        xb_0=xb_0+dxb;
  }
}

void funktion_b_ohne(double xb_0, double xb_f, int n, double* ybo){ // Funktion aus Aufgabenteil a
    double dxb = (xb_f-xb_0)/n;
    for(int i=0; i<=n; i++){
        ybo[i] = (sin(xb_0))/(1+cos(xb_0));
        xb_0=xb_0+dxb;
  }
}

void funktion_c(double xc_0, double xc_f, double delta, int n, double* yc){ // Funktion aus Aufgabenteil a
    double dxc = (xc_f-xc_0)/n;
    for(int i=0; i<=n; i++){
        yc[i] = sin(xc_0+delta)-sin(xc_0);
        xc_0=xc_0+dxc;
  }
}



int main(){

  cout << "Start Aufgabe 2" << "\n";

  double xa_0 = 10000000000000.0; // Festlegung des Startwerts Aufgabenteil a
  double xa_f = 1000000000000000.0 ; // Festlegung der Endwertes Aufgabenteil a

  double xb_0 = 0.000000001; // Festlegung des Startwerts Aufgabenteil a
  double xb_f = 0.0000001 ; // Festlegung der Endwertes Aufgabenteil a

  double xc_0 = 0.0; // Festlegung des Startwerts Aufgabenteil a
  double xc_f = 3.0 ; // Festlegung der Endwertes Aufgabenteil a

  double delta= 0.0000001;

  int n = 100.0; // Anzahl der Schritte

  cout << "Anzahl der Schritte beträgt " << n << "\n";

  double y_a[n];
  funktion_a(xa_0, xa_f, n, y_a); // Aufruf des normalen Eulerverfahrens

  double y_b[n];
  funktion_b(xb_0, xb_f, n, y_b); // Aufruf des normalen Eulerverfahrens

  double y_bo[n];
  funktion_b_ohne(xb_0, xb_f, n, y_bo); // Aufruf des normalen Eulerverfahrens

  double y_c[n];
  funktion_c(xc_0, xc_f, delta, n, y_c); // Aufruf des normalen Eulerverfahrens



  ofstream myfile_a;
  myfile_a.open("./Data/two_a.txt");
  myfile_a << "# x_a f_a(xa) f_a_ohne(xa)" << "\n";
  for(int i=0; i<=n; i++){
    myfile_a << xa_0+i*(xa_f-xa_0)/n << " ";
    myfile_a << y_a[i] << "\n";
    //myfile << y_ao[i] << " ";
  }
  myfile_a.close();



  ofstream myfile_b;
  myfile_b.open("./Data/two_b.txt");
  myfile_b << "# x_b f_b(xb) f_b_ohne(xb)" << "\n";
  for(int i=0; i<=n; i++){
    myfile_b << xb_0+i*(xb_f-xb_0)/n << " ";
    myfile_b << y_b[i] <<  "\n";
    myfile_b << y_bo[i] <<  "\n";
  }
  myfile_b.close();

cout << "Hier bin ich " << n << "\n";

  ofstream myfile_c;
  myfile_c.open("./Data/two_c.txt");
  myfile_c << "# x_c f_c(xb) f_c_ohne(xc)" << "\n";
  for(int i=0; i<=n; i++){
    myfile_c << xc_0+i*(xc_f-xc_0)/n << " ";
    myfile_c << y_c[i] << "\n";
    //myfile << y_co[i] << " ";
  }
  myfile_c.close();

  cout << "Aufgabe 2 ist abgeschlossen" << "\n";
  return 0;
}
