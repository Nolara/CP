#include <iostream>
#include <math.h> // exp-Funktion
#include <fstream> // Ausgabe als txt
#include <vector>
using namespace std;

vector<double> euler_normal(int n, double y_0, double delta_t){ //Definiere normales Eulerverfahren
    vector<double> y;
    y.push_back(y_0);
    for(int i=1; i<=n; i++){
        y.push_back(y[i-1]*(1-delta_t));
    }
    return y;
}

vector<double> euler_symm(int n, double y_0,  double y_1, double delta_t){ //Definiere symmetrisches Eulerverfahren
    vector<double> y;
    y.push_back(y_0);
    y.push_back(y_1);
    for(int i=2; i<=n; i++){
        y.push_back(-2*delta_t*y[i-1]+y[i-2]);
    }
    return y;
}

int main(){

  cout << "Start Aufgabe 3" << endl;

  double delta_t = 0.1  ; // Festlegung der Schrittweite
  double y_0a = 1.0; // Festlegung des Anfangswerts Aufgabenteil a
  double y_0b = 1.0-delta_t; // Festlegung des Anfangswerts Aufgabenteil b
  double y_1a = exp(-delta_t); // Festlegung des Anfangswerts Aufgabenteil a
  double y_1b = y_0a-delta_t; // Festlegung des Anfangswerts Aufgabenteil b
  int n = 10/delta_t; // Anzahl der Schritte

  vector<double> y_a_normal;
  y_a_normal = euler_normal(n, y_0a, delta_t); // Aufruf des normalen Eulerverfahrens Aufgabenteil a

  vector<double> y_a_symm;
  y_a_symm = euler_symm(n, y_0a, y_1a, delta_t); // Aufruf des symmetrischen Eulerverfahrens Aufgabenteil a

  vector<double> y_b_normal;
  y_b_normal = euler_normal(n, y_0b, delta_t); // Aufruf des normalen Eulerverfahrens Aufgabenteil b

  vector<double> y_b_symm;
  y_b_symm = euler_symm(n, y_0a, y_1b, delta_t); // Aufruf des symmetrischen Eulerverfahrens Aufgabenteil b


  ofstream afile ("Data/three_a.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil a
  afile << "# t=n*delta_t euler_normal euler_symmetrisch analytisch" << endl;
  for(int i=0; i<=n; i++){
    // Schreibe Daten in txt Datei
    afile << i*delta_t << " ";
    afile << y_a_normal[i] << " ";
    afile << y_a_symm[i] << " ";
    afile << exp(-i*delta_t) << endl; // Analytische Lösung
  }
  afile.close();

  ofstream bfile ("Data/three_b.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil b
  bfile << "# t=n*delta_t euler_normal euler_symmetrisch analytisch" << endl;
  for(int i=0; i<=n; i++){
    // Schreibe Daten in txt Datei
    bfile << i*delta_t << " ";
    bfile << y_b_normal[i] << " ";
    bfile << y_b_symm[i] << " ";
    bfile << exp(-i*delta_t) << endl; // Analytische Lösung
  }
  bfile.close();

  cout << "Aufgabe 3 ist abgeschlossen" << endl;
  return 0;
}
