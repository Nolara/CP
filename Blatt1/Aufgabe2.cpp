#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{


   VectorXf x(10);
   x <<0, 2.5, -6.3, 4, -3.2, 5.3, 10.1, 9.5, -5.4, 12.7; //x aus der Aufgabenstellung

   VectorXf a(10);
   a <<1, 1, 1, 1, 1, 1, 1, 1, 1, 1; //Hilfsvektor a

   MatrixXf A(10,2);
   A << x, a; //Erstelle Matrix aus x und a

   VectorXf y(10);
   y <<4, 4.3, -3.9, 6.5, 0.7, 8.6, 13, 9.9, -3.6, 15.1; //y aus der Aufgabenstellung

   MatrixXf AT(2,10);
   AT << A.transpose(); //Transponiere A

   MatrixXf ATA(2,2);
   ATA << AT*A; //Produkt aus A transponiert und A

   VectorXf ATy(2);
   ATy << AT*y; //Produkt aus A transponiert und y

   Vector2f x0;
   x0= ATA.partialPivLu().solve(ATy); //Löse ATAx= ATy

   ofstream afile ("Data/two_data.txt", std::ofstream::out); // Erstelle txt Datei für Daten aus Aufgabenstellung
   ofstream bfile ("Data/two_x0.txt", std::ofstream::out); // Erstelle txt Datei Lösungsvektor
   afile << "# x y" << "\n";
   bfile << "# m n" << "\n";
   int n;
   n=x.size();
   for(int i=0; i<n; i++){
       afile << x[i] << " ";
       afile << y[i] << "\n";
   }
   bfile << x0[0] << " " << x0[1] << "\n";

   afile.close();
   bfile.close();

}
