#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{

   Vector3f a1;
   Vector3f a2;
   Vector3f a3;

   a1 <<0.5, sqrt(3)/2, 0;
   a2 <<-0.5, sqrt(3)/2, 0;
   a3 <<0, 0, 1;

   //A << 0.5,-0.5,0,  sqrt(3)/2,sqrt(3)/2,0,  0,0,1;

   Vector3f b;
   b << 2, 0, 2;

   Matrix3f A;
   A << a1, a2, a3;

   Vector3f x;
   x= A.partialPivLu().solve(b);

   Matrix3f P;
   P = A.partialPivLu().permutationP();

   Matrix3f LU;
   LU=A.partialPivLu().matrixLU();

   Matrix3f L;
   L=LU.triangularView<StrictlyLower>();
   for(int i=0; i<=2; i++){
     L(i,i)=1;
   }

   Matrix3f U;
   U=LU.triangularView<Upper>();

   Vector3f b2;
   b2 << 1, 2*sqrt(3), 3;

   Vector3f pb;
   pb=P*b2;

   Vector3f x_zwischen;

   x_zwischen(0)=pb(0)/L(0,0);

   for(int i = 1; i < 3; i++)
   {
       float s = 0;
       for(int j = 1; j < i; j++)
       {
           s = s + L(i,j)*x_zwischen(j);
       }
       x_zwischen(i) = (pb(i) - s)/L(i,i);
    }

   Vector3f x2;

   x2(2)=x_zwischen(2)/U(2,2);

   for (int i = 1; i >=0; i--)
   {
       float s = 0;
       for (int j = i+1; j <3; ++j)
       {
           s = s + U(i,j)*x2(j);;
       }
       x2(i)=(x_zwischen(i)-s)/U(i,i);
   }

   Matrix3f A2;
   A2 << a3, a2, a1;

   Vector3f xd;
   xd= A2.partialPivLu().solve(b);

   Matrix3f P2;
   P2 = A2.partialPivLu().permutationP();

   Matrix3f LU2;
   LU2=A2.partialPivLu().matrixLU();

   Matrix3f L2;
   L2=LU2.triangularView<StrictlyLower>();
   for(int i=0; i<=2; i++){
     L2(i,i)=1;
   }

   Matrix3f U2;
   U2=LU2.triangularView<Upper>();

   ofstream afile ("Data/one.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil a
   afile << "# Aufgabenteil b:" << "\n" <<  "\n";
   afile << "# Lösungsvektor:" << "\n";
   int n_b;
   n_b=x.size();
   for(int i=0; i<n_b; i++){
       afile << x[i] << "\n";
   }
   afile << "\n" << "# Matrix P:" << "\n";

   int r=P.rows();;
   int c=P.cols();

   for (int i = 0; i < r; ++i)
   {
       for (int j = 0; j < c; ++j)
       {
           afile << P(i,j) << " ";
       }
       afile << "\n";
   }

    afile << "\n" << "# Matrix L:" << "\n";

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << L(i,j) << " ";
        }
        afile << "\n";
    }

    afile << "\n" << "# Matrix U:" << "\n";

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << U(i,j) << " ";
        }
        afile << "\n";
    }

    afile << "\n" << "# Aufgabenteil c:" << "\n" << "\n";
    afile << "# Lösungsvektor:" << "\n";

    for(int i=0; i<n_b; i++){
        afile << x2[i] << "\n";
    }

    afile << "\n" <<  "# Aufgabenteil d:" << "\n" << "\n";
    afile << "# Lösungsvektor:" << "\n";

    for(int i=0; i<n_b; i++){
        afile << xd[i] << "\n";
    }
    afile << "\n" << "# Matrix P:" << "\n";

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            afile << P2(i,j) << " ";
        }
        afile << "\n";
    }

     afile << "\n" << "# Matrix L2:" << "\n";

     for (int i = 0; i < r; ++i)
     {
         for (int j = 0; j < c; ++j)
         {
             afile << L(i,j) << " ";
         }
         afile << "\n";
     }

     afile << "\n" << "# Matrix U2:" << "\n";

     for (int i = 0; i < r; ++i)
     {
         for (int j = 0; j < c; ++j)
         {
             afile << U(i,j) << " ";
         }
         afile << "\n";
     }


   afile.close();

}
