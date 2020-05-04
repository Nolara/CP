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

//Funktion zum Hinzufügen einer Zeile mit Inhalten "content" zur txt Datei mit Namen "name"
void addcontent(const std::string& name, long long& content, bool append = false) {
    std::ofstream outfile;
    if (append)
        outfile.open(name, std::ios_base::app);
    else
        outfile.open(name);
    outfile << content << "\n";
}

//Funktion zur Messung der Zeit der Invertierung mit Rückgabe des Ergebnissvektors
VectorXd Invert(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 0 );
    MatrixXd Inv(size,size);
    Inv=M.inverse(); //Matrix Invertieren
    x= Inv*b; //Gleichungssystem lösen
    Profiler::stop( 0 );
    long long time_inv = Profiler::getTimeInNS( 0 );
    addcontent("Data/three_t_inverse.txt", time_inv ,true); //Wert des Timers zur txt Datei hinzufügen
    return x;
}

//Funktion zur Messung der Zeit der LU Zerlegung  mit vollständiger Pivotisierung mit Rückgabe des Ergebnissvektors
VectorXd FullLu(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 1 );
    MatrixXd FullLu(size,size);
    FullPivLU<MatrixXd> flu(M); //LU Zerlegung
    x= flu.solve(b); //Lösen des LGS
    Profiler::stop( 1 );
    long long time_full = Profiler::getTimeInNS( 1 );
    addcontent("Data/three_t_full_LU.txt", time_full ,true); //Wert des Timers zur txt Datei hinzufügen

    //Errechne Rang der Matrix
    int r;
    r= flu.rank();

    //Prüfe ob Matrix invertierbar
    if (r!=size)
    {
        cout << "Matrix ist nicht invertierbar!";
    }

    return x;
}

//Funktion zur Messung der Zeit der LU Zerlegung  mit teilweiser Pivotisierung mit Rückgabe des Ergebnissvektors
VectorXd PartLu(int size, MatrixXd M, VectorXd b){
    VectorXd x(size);
    Profiler::start( 2 );
    MatrixXd PartLu(size,size);
    PartialPivLU<MatrixXd> plu(M); //LU Zerlegung
    x= plu.solve(b); //Lösen des LGS
    Profiler::stop( 2 );
    long long time_part = Profiler::getTimeInNS( 2 );
    addcontent("Data/three_t_partial_LU.txt", time_part ,true); //Wert des Timers zur txt Datei hinzufügen
    return x;
}


int main()
{
    int steps =13; //Festlegung der Anzahl der Schriite

    //txt Dateien erstellen
    ofstream afile ("Data/three_t_inverse.txt", std::ofstream::out);
    ofstream bfile ("Data/three_t_full_LU.txt", std::ofstream::out);
    ofstream cfile ("Data/three_t_partial_LU.txt", std::ofstream::out);
    afile.close();
    bfile.close();
    cfile.close();

    Profiler::init( 3 ); //Erstelle drei Timer

    //Erstelle Matrix für relativen Fehler der Lösungen
    MatrixXd Abw(steps,3);

    for (int k = 1; k < steps; ++k)
    {
        //Verwende logarithmische Schritte
        int i;
        i=pow(2,k);
        //Erstelle Zufallsmatrix und Zufallsvektor
        MatrixXd M(i,i);
        VectorXd b(i);
        M=Random_matrix(i);
        b=Random_vector(i);

        //Rufe Funktion zur Berechnung mit vollständig pivotisierter LU Zerlegung auf
        VectorXd x_FullLU(i);
        x_FullLU = FullLu(i, M, b);

        //Rufe Funktion zur Berechnung mit Inverser auf
        VectorXd x_Inv(i);
        x_Inv = Invert(i, M, b);

        //Rufe Funktion zur Berechnung mit teilweise pivotisierter LU Zerlegung auf
        VectorXd x_parLU(i);
        x_parLU = PartLu(i, M, b);

        //Berechne relative Abweichung
        Abw(k,0)=(M*x_Inv - b).norm() / b.norm();
        Abw(k,1)=(M*x_FullLU - b).norm() / b.norm();
        Abw(k,2)=(M*x_parLU - b).norm() / b.norm();

        Profiler::resetAll(); //Setze alle Timer zurück
    }
    //Speichere relative Abweichung in txt Datei
    ofstream dfile ("Data/three_relative_error.txt", std::ofstream::out);
    dfile << "# relativer Fehler Invertierung, relativer Fehler Full LU, relativer Fehler Partial LU " << "\n" ;
    dfile << Abw;
    dfile.close();
}
