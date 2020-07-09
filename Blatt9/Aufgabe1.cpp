#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include <random>

using namespace std;
using namespace Eigen;


//Metropolis Algorithmus
double Metropolis(double N_steps, double Magnetfeld, mt19937 &generator, uniform_real_distribution<double> &distribution ){

  VectorXd spin = VectorXd::Zero(N_steps); //Spinvektor für +1/-1
  double m=0.0; //Summe über Spinausrichtungen = Magnetisierung

  //belieber Startzustand i_0:
  uniform_int_distribution<int> startdistribution(0, 1); // liefert den Wert 0 oder 1
  spin[0] = startdistribution(generator)*2 - 1; // liefert den Wert +1 oder -1 für den Spin

  for(int n=0; n<N_steps-1;n++){

      //MC vorschlagen:
      double p = distribution(generator);
      double delta_E = 2*spin[n]*Magnetfeld;
      bool accept;

			//Überprüfen, ob der MC Move akzepiert wird
      if(delta_E<0){
        accept = true;
      }
      else{
        if(p<exp(-delta_E)){
          accept = true;
        }
        else{
          accept = false;
        }
      }

      //k->k+1, falls move akzepiert ist
      if(accept==true){
        spin[n+1]=(-1)*spin[n];
      }
      else{
        spin[n+1]=spin[n];
      }

  }
  m=spin.sum()/N_steps;
  return m;
}


int main(){

  random_device zufall;
  double N =1e5;
  mt19937 generator(zufall());
  uniform_real_distribution<double> distribution(0,1);

  VectorXd Magnetfeld = VectorXd::LinSpaced(1e4, -5, 5);
  ofstream file;
  string filename = "Data/one.txt";
  file.open(filename.c_str());
  file << "# Magnetfeld, Magnetisierung m \n";
  for(int i=0; i<1e4; i++){
    file << Magnetfeld[i] << "\t" << Metropolis(N,Magnetfeld[i],generator,distribution) << "\n" ;
  }
  file.close();

  return 0;
}
