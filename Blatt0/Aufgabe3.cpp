#include <iostream>
#include <math.h> // exp-Funktion
#include <fstream> // Ausgabe als txt
using namespace std;


void euler_normal(int n, double y_0, double delta_t, double* y){ // Normales Eulerverfahren
  y[0] = y_0; //Startwertübergabe
  for(int i=1; i<=n; i++){
    y[i] = y[i-1]*(1-delta_t); //Rekursion
  }
}


void euler_symm(int n, double y_0, double delta_t, double* y){
  y[0] = y_0; // Startwertübergabe
  y[1] = exp(-delta_t);
  for(int i=2; i<=n; i++){
    y[i] = -2*delta_t*y[i-1]+y[i-2]; //Rekursion
  }
}


int main(){

  cout << "Start Aufgabe 3" << endl;
  double y_0 = 1.0; // Festlegung des Startwerts
  double delta_t = 0.1; // Festlegung der Schrittweite
  int n = 10/delta_t; // Anzahl der Schritte
  cout << "Anzahl der Schritte beträgt " << n << endl;
  double y_normal[n]; //Werte des normalen Eulerverfahrens
  euler_normal(n, y_0, delta_t, y_normal); // Aufruf des normalen Eulerverfahrens
  double y_symm[n]; //Werte des symmetrischen Eulerverfahrens
  euler_symm(n, y_0, delta_t, y_symm); // Aufruf des symmetrischen Eulerverfahrens

  ofstream myfile;
  myfile.open("./Data/three.txt");
  myfile << "# t=n*delta_t euler_normal euler_symmetrisch analytisch" << endl;
  for(int i=0; i<=n; i++){
    myfile << i*delta_t << " ";
    myfile << y_normal[i] << " ";
    myfile << y_symm[i] << " ";
    myfile << exp(-i*delta_t) << endl; // Analytische Lösung
  }
  myfile.close();

  cout << "Aufgabe 3 ist abgeschlossen" << endl;
  return 0;
}
