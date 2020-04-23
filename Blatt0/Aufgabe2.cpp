#include <iostream>
#include <math.h> // exp-Funktion
#include <fstream> // Ausgabe als txt
#include <vector>
using namespace std;


vector<float> funktion_a(float xa_0, float xa_f, int n){ // Funktion aus Aufgabenteil a
    float dxa = (xa_f-xa_0)/n;
    vector<float> ya;
    for(int i=0; i<=n; i++){
        ya.push_back(1/sqrt(xa_0)-1/sqrt(xa_0+1));
        xa_0=xa_0+dxa;
    }
    return ya;
}

vector<float> better_a(float xa_0, float xa_f, int n){ // Funktion aus Aufgabenteil a
    float dxa = (xa_f-xa_0)/n;
    vector<float> za;
    for(int i=0; i<=n; i++){
        za.push_back(1/((xa_0+1)*sqrt(xa_0)+xa_0*sqrt(xa_0+1)));
        xa_0=xa_0+dxa;
    }
    return za;
}

vector<float> funktion_b(float xb_0, float xb_f, int n){
    float dxb = (xb_f-xb_0)/n;
    vector<float> yb;
    for(int i=0; i<=n; i++){
        yb.push_back((1-cos(xb_0))/sin(xb_0));
        xb_0=xb_0+dxb;
    }
    return yb;
}

vector<float> better_b(float xb_0, float xb_f, int n){
    float dxb = (xb_f-xb_0)/n;
    vector<float> zb;
    for(int i=0; i<=n; i++){
        zb.push_back(sin(xb_0)/(1+cos(xb_0)));
        xb_0=xb_0+dxb;
    }
    return zb;
}

vector<float> funktion_c(float xc_0, float xc_f, float delta, int n){ // Funktion aus Aufgabenteil a
    float dxc = (xc_f-xc_0)/n;
    vector<float> yc;
    for(int i=0; i<=n; i++){
        yc.push_back(sin(xc_0+delta)-sin(xc_0));
        xc_0=xc_0+dxc;
    }
    return yc;
}

vector<float> better_c(float xc_0, float xc_f, float delta, int n){ // Funktion aus Aufgabenteil a
    float dxc = (xc_f-xc_0)/n;
    vector<float> zc;
    for(int i=0; i<=n; i++){
        zc.push_back(cos(xc_0)*sin(delta)-(sin(xc_0)*sin(delta)*sin(delta))/(cos(delta)+1));
        xc_0=xc_0+dxc;
    }
    return zc;
}



int main(){

  cout << "Start Aufgabe 2" << "\n";

  float xa_0 = 100000.0; // Festlegung des Startwerts Aufgabenteil a
  float xa_f = 10000000.0 ; // Festlegung der Endwertes Aufgabenteil a

  float xb_0 = 0.00001; // Festlegung des Startwerts Aufgabenteil a
  float xb_f = 0.002 ; // Festlegung der Endwertes Aufgabenteil a

  float xc_0 = 0.0; // Festlegung des Startwerts Aufgabenteil a
  float xc_f = 4.0 ; // Festlegung der Endwertes Aufgabenteil a

  float delta= 0.000001;

  int n = 10000.0; // Anzahl der Schritte

  cout << "Anzahl der Schritte beträgt " << n << "\n";

  vector<float> y_a;
  y_a = funktion_a(xa_0, xa_f, n); // Aufruf des normalen Eulerverfahrens

  vector<float> z_a;
  z_a = better_a(xa_0, xa_f, n); // Aufruf des normalen Eulerverfahrens

  vector<float> y_b;
  y_b =funktion_b(xb_0, xb_f, n); // Aufruf des normalen Eulerverfahrens

  vector<float> z_b;
  z_b = better_b(xb_0, xb_f, n); // Aufruf des normalen Eulerverfahrens

  vector<float> y_c;
  y_c = funktion_c(xc_0, xc_f, delta, n); // Aufruf des normalen Eulerverfahrens

  vector<float> z_c;
  z_c = better_c(xc_0, xc_f, delta, n); // Aufruf des normalen Eulerverfahrens


  ofstream afile ("Data/two_a.txt", std::ofstream::out);
  ofstream bfile ("Data/two_b.txt", std::ofstream::out);
  ofstream cfile ("Data/two_c.txt", std::ofstream::out);
  ofstream del ("Data/two_delta.txt", std::ofstream::out);
  afile << "# x_a f_a(xa) f_a_ohne(xa)" << "\n";
  bfile << "# x_b f_b(xb) f_b_ohne(xb)" << "\n";
  cfile << "# x_c f_c(xc) f_c_ohne(xc)" << "\n";
  del << "# delta " << "\n";
  del << delta << "\n";
  for(int i=0; i<=n; i++){
      afile << xa_0+i*(xa_f-xa_0)/n << " ";
      bfile << xb_0+i*(xb_f-xb_0)/n << " ";
      cfile << xc_0+i*(xc_f-xc_0)/n << " ";
      afile << y_a[i] << " ";
      bfile << y_b[i] << " ";
      cfile << y_c[i] << " ";
      afile << z_a[i] << "\n";
      bfile << z_b[i] << "\n";
      cfile << z_c[i] << "\n";
  }

  afile.close();
  bfile.close();
  cfile.close();

  cout << "Aufgabe 2 ist abgeschlossen" << "\n";
  return 0;
}
