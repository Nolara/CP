#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include "profiler.h"

using namespace std;
using namespace Eigen;

//Funktion zum Erstellen einer Zufallsmatrix der Größe nxn
MatrixXd Random_matrix(int n){
    MatrixXd A(n,n);
    A= MatrixXd::Random(n,n);
    return A;
}

//Funktion zum Erstellen eines Zufallsvektors der Größe n
VectorXd Random_vector(int n){
    VectorXd v(n);
    v=VectorXd::Random(n);
    return v;
}

//Funktion zur Messung der Zeit der einzelnen Schritte bis zur maximalen Dimension steps
MatrixXd Zeit(int steps){
    MatrixXd Zeit(steps, 4);
    Profiler::init( 4 ); //Erstelle 4 Timer, je einer für jeden Schritt
    for (int k = 0; k < steps; ++k)
    {
        //Verwende logarithmische Schritte
        int i;
        i=pow(2,k);

        //Messung der Zeit zum Erstellen der Zufallsmatrix
        Profiler::start( 0 );
        MatrixXd A=Random_matrix(i);
        Profiler::stop( 0 );
        Zeit(k,0)=Profiler::getTimeInNS( 0 );

        //Messung der Zeit zum Erstellen des Zufallsvektors
        Profiler::start( 1 );
        VectorXd b=Random_vector(i);
        Profiler::stop( 1 );
        Zeit(k,1)=Profiler::getTimeInNS( 1 );

        //Messung der Zeit der LU Zerlegung
        Profiler::start( 2 );
        PartialPivLU<MatrixXd> lu(A);
        Profiler::stop( 2 );
        Zeit(k,2)=Profiler::getTimeInNS( 2 );

        //Messung der Zeit zum Lösen des Gleichungssystems
        Profiler::start( 3 );
        VectorXd x(i);
        x= lu.solve(b);
        Profiler::stop( 3 );
        Zeit(k,3)=Profiler::getTimeInNS( 3 );

        Profiler::resetAll(); //Alle Timer zurücksetzen
    }
    return Zeit;
}


int main()
{
    int steps =13; //Festlegung der Anzahl der schritte
    MatrixXd Results(steps,4);
    Results = Zeit(steps); //Aufrufen der Funktion Zeit

    //Abspeichern der Ergebnisse
    ofstream afile ("Data/two_data.txt", std::ofstream::out);
    afile << "# Erzeugung Zufallsmatrix, Erzeugung Zufallsvektor, LU Zerlegung, Lösung" << "\n";
    afile << Results;
    afile.close();

}
