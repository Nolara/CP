#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{

   Vector3d a1;
   Vector3d a2;
   Vector3d a3;

   //Eingabe der Vektoren aus der Aufgabenstellung
   a1 <<0.5, sqrt(3)/2, 0;
   a2 <<-0.5, sqrt(3)/2, 0;
   a3 <<0, 0, 1;

   Vector3d b;
   b << 2, 0, 2;

   //Definiere Matrix aus Basisvektoren
   Matrix3d A;
   A << a1, a2, a3;

   //Lösen mit LU Zerlegung
   Vector3d x;
   x= A.partialPivLu().solve(b);

   //Bestimme Matrix Pivotisierungsmatrix P
   Matrix3d P;
   P = A.partialPivLu().permutationP();

   //Bestimme Matrix L*U
   Matrix3d LU;
   LU=A.partialPivLu().matrixLU();

   //Teile L*U in L und U auf
   Matrix3d L;
   L=LU.triangularView<StrictlyLower>();
   for(int i=0; i<=2; i++){ //Setze Diagonale =1
     L(i,i)=1;
   }

   Matrix3d U;
   U=LU.triangularView<Upper>();

   //Eingabe des Vektors aus Aufgabenteil c)
   Vector3d b2;
   b2 << 1, 2*sqrt(3), 3;

   //Multipliziere b mit Pivotisierungsmatrix P
   Vector3d pb;
   pb=P*b2;

   //Erstelle Vektor zum Speichern der Zwischenergebnisse aus der Forward Substitution
   Vector3d x_zwischen;

   //Erster Wert der Forward Substitution
   x_zwischen(0)=pb(0)/L(0,0);

   //Forward Substitution
   for(int i = 1; i < 3; i++)
   {
       float s = 0;
       for(int j = 0; j < i; j++)
       {
           s = s + L(i,j)*x_zwischen(j);
       }
       x_zwischen(i) = (pb(i) - s)/L(i,i);
   }

   //Definiere Lösungsvektor für c)
   Vector3d x2;

   //Erster Wert der Backward Substitution
   x2(2)=x_zwischen(2)/U(2,2);

   //Backward Substitution
   for (int i = 1; i >=0; i--)
   {
       double s = 0;
       for (int j = i+1; j <3; ++j)
       {
           s = s + U(i,j)*x2(j);;
       }
       x2(i)=(x_zwischen(i)-s)/U(i,i);
   }


   //Umgestellte Matrix aus Aufgabenteil d)
   Matrix3d A2;
   A2 << a3, a2, a1;

   //Löse d) mit LU
   Vector3d xd;
   xd= A2.partialPivLu().solve(b);

   //Pivotisierungsmatrix P aus d)
   Matrix3d P2;
   P2 = A2.partialPivLu().permutationP();

   //Matrix LU aus d)
   Matrix3d LU2;
   LU2=A2.partialPivLu().matrixLU();

   //Matrix L aus d)
   Matrix3d L2;
   L2=LU2.triangularView<StrictlyLower>();
   for(int i=0; i<=2; i++)
   {
       L2(i,i)=1;
   }

   //Matrix U aus d)
   Matrix3d U2;
   U2=LU2.triangularView<Upper>();

   //Schreibe Werte in txt Datei
   ofstream afile ("Data/one.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil a
   afile << "# Aufgabenteil b:" << "\n" <<  "\n";
   afile << "# Lösungsvektor:" << "\n";
   int n_b;
   n_b=x.size(); //Länge Vektor x
   for(int i=0; i<n_b; i++){
       afile << x[i] << "\n"; //Schreibe x in txt Datei
   }
   afile << "\n" << "# Matrix P:" << "\n";

   //Größe der Matrizen
   int r=P.rows();;
   int c=P.cols();

   //Schreibe P in txt Datei
   for (int i = 0; i < r; ++i)
   {
       for (int j = 0; j < c; ++j)
       {
           afile << P(i,j) << "\t";
       }
       afile << "\n";
   }

   //Schreibe L in txt Datei
    afile << "\n" << "# Matrix L:" << "\n";

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << L(i,j) << "\t";
        }
        afile << "\n";
    }

    //Schreibe U in txt Datei
    afile << "\n" << "# Matrix U:" << "\n";

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << U(i,j) << "\t";
        }
        afile << "\n";
    }

    afile << "\n" << "# Aufgabenteil c:" << "\n" << "\n";
    afile << "# Lösungsvektor:" << "\n";

    //Schreibe x aus c) in txt Datei
    for(int i=0; i<n_b; i++){
        afile << x2[i] << "\n";
    }

    afile << "\n" <<  "# Aufgabenteil d:" << "\n" << "\n";
    afile << "# Lösungsvektor:" << "\n";

    //Schreibe x aus d) in txt Datei
    for(int i=0; i<n_b; i++){
        afile << xd[i] << "\n";
    }
    afile << "\n" << "# Matrix P:" << "\n";

    //Schreibe P aus d) in txt Datei
    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << P2(i,j) << "\t";
        }
        afile << "\n";
    }

    //Schreibe L aus d) in txt Datei
     afile << "\n" << "# Matrix L2:" << "\n";

     for (int i = 0; i < r; ++i)
     {
         for (int j = 0; j < c; ++j)
         {
             afile << L2(i,j) << "\t";
         }
         afile << "\n";
     }

     //Schreibe U aus d) in txt Datei
     afile << "\n" << "# Matrix U2:" << "\n";

     for (int i = 0; i < r; ++i)
     {
         for (int j = 0; j < c; ++j)
         {
             afile << U2(i,j) << "\t";
         }
         afile << "\n";
     }


   afile.close();

}
