#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "service.cpp"

using namespace std;
using namespace Eigen;



//Funktion zum Erstellen des Namens der txt Datzei
string name( const string& basename, float index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

//Funktion zur Erstellung der Rang-k Reduktion
void Reduktion(int k, MatrixXd Singular, MatrixXd U, MatrixXd V, int n)
{
    MatrixXd Reduktion(n, n);

    for (int i = 0; i < k; i++)
    {
        //Die ersten k Singulärwerte werden zusammen mit den entsprechenden Eigenvektoren
        //zur Matrix "Reduktion hinzugefügt"
    	Reduktion += Singular(i) * U.col(i) * V.row(i);
    }
    //Erstelle txt Datei für ebtsprechendes k und speichere Matrix ab
    ofstream file;
    file.open (name("Data/One_Reduktion_", k, ".txt"), ofstream::out);
    file << Reduktion.transpose();
    file.close();
    Reduktion.resize(0,0);
}

int main()
{
    int size = 512; //Größe der Bilddatei
    MatrixXd M(size, size);
    loadData(M, "Bild", size, size); //Bilddatei einlesen
    MatrixXd Original(size, size);
    Original = M.transpose(); //Originalbild zum Vergleich abspeichern

    //Singulärwertzerlegung berechnen
    BDCSVD<MatrixXd> svd(M, ComputeFullV | ComputeFullU );

    //Singulärwerte in Matrix abspeichern
    MatrixXd Singular(size, size);
    Singular=svd.singularValues();

    //Matrix U
    MatrixXd U(size, size);
    U=svd.matrixU();

    //Matrix V
    MatrixXd V(size, size);
    V=svd.matrixV().transpose();

    //Originalbild abspeichern
    ofstream afile ("Data/one_original.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil a
    afile << "# Aufgabenteil b:" << "\n" <<  "\n";
    afile << Original << endl;
    afile.close();

    //Funktion "Reduktion" für verschiedene k aufrufen
    Reduktion(10, Singular,U, V, size);
    Reduktion(20, Singular,U, V, size);
    Reduktion(50, Singular,U, V, size);

}
