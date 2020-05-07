#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include "profiler.h"

using namespace std;
using namespace Eigen;

int main()
{
    int Schritte = 14;
    MatrixXd Zeit(Schritte,4);
    Profiler::init( 3 );

    for(int N=0; N<Schritte; N++)
    {
      int logN = pow(2,N);
      VectorXd b(logN);
      b=VectorXd::Random(logN);         

      Profiler::start( 0 );        //Laufzeitmessung für random NxN Matrix 
      MatrixXd M(logN,logN);
      M =MatrixXd::Random(logN,logN);       
      Profiler::stop( 0 );

      Profiler::start( 1 );         //Laufzeitmessung LU-Zerlegung mit partialer Pivotierung
      PartialPivLU<MatrixXd> LU(M);
      Profiler::stop(1);

      Profiler::start( 2 );         //Laufzeitmessung für Lösung des Gleichungssystems
      VectorXd x(logN);
      x = LU.solve(b);
      Profiler::stop( 2 );

      Zeit(N,0) = logN;
      Zeit(N,1) = Profiler::getTimeInNS( 0 );
      Zeit(N,2) = Profiler::getTimeInNS( 1 );
      Zeit(N,3) = Profiler::getTimeInNS( 2 );
      Profiler::resetAll();
    }

    ofstream filea ("two.txt", std::ofstream::out);    // Speichern der Laufzeiten in Abhängigkeit von N
    filea << "#N   Laufzeiten für 1)   2)    3)" << '\n' ;
    filea << Zeit ;
    filea.flush();
    filea.close();

    return 0;
}
