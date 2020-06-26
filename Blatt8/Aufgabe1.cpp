#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>


using namespace std;
using namespace Eigen;

//Funktion zur Namensgebung der erstellten txt Dateien
string give_name( const string& basename, string index, const string& ext) {
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

//Kraft durch das Lennard Jones Potential
//Übergabewerte sind der Vektor r zwischen 2 Teilchen und der kritische Abstand rc
//Es wird der Kraftvektor ausgegeben
VectorXd Kraft_LJ(VectorXd &r, double &rc){
	const int dim=r.size();
	double r_squared=r.squaredNorm();


	if (r_squared>(rc*rc)) {
		VectorXd f=VectorXd::Zero(dim);
		return f;
	} else if (r_squared ==0) {
		r_squared+=1e-9;
		VectorXd f=24*(2*pow(r_squared,-7)-pow(r_squared,-4))*r;
		return f;
	} else {
		VectorXd f=24*(2*pow(r_squared,-7)-pow(r_squared,-4))*r;
		return f;
	}

}

//Lennard Jones Potential zur Berechnung der Potentiellen Energie
//Übergabewerte sind der Vektor r zwischen 2 Teilchen und der kritische Abstand rc
//Es wird der Werte der potentiellen Energie ausgegeben
double Pot_LJ(VectorXd &r, double &rc){
	double r_squared=r.squaredNorm();

	if (r_squared>(rc*rc)) {
		double f=0;
		return f;
	} else if (r_squared ==0) {
		r_squared+=1e-9;
		double f=4*(pow(r_squared,-6)-pow(r_squared,-3));
		return f;
	} else {
		double f=4*(pow(r_squared,-6)-pow(r_squared,-3));
		return f;
	}


}

//Funktion zur Berechnung der Beschleunigung
//Übergabewerte sind:  die Funktion f, die die Kraft beschreibt, die Matrix Yn mit den Orts und Geschwindigkeitsvektoren,
//die Größe L des Kastens, die Funktion pot zur Berechnung der potentiellen Energie und ein double Epot, in dem die
// potentielle Energie gespeichert wird
//Es wird die Matrix mit der Beschleunigung ausgegeben
MatrixXd accelaration(VectorXd (*f)(VectorXd &, double &), MatrixXd &Yn, double L, double (*pot)(VectorXd &, double &), double &Epot) {
	int row= Yn.rows();
	int col= Yn.cols();
	double rc=L/2.0;
	//Verschiebungen aller Teilchen um +-L in beiden Dimensionen um virtuelle Teilchen zu simulieren
	MatrixXd A=MatrixXd::Zero(row/2,col);

	for (int x = 0; x < col; ++x)
	{
		for (int xstrich = 0; xstrich < col; ++xstrich)
		{
			if (x!=xstrich)
			{
				VectorXd Verschiebung=VectorXd::Zero(row/2);
				for (int i = -1; i <= 1; ++i)
				{
					Verschiebung(0)=i*L;
					for (int j = -1; j <= 1; ++j)
					{
						Verschiebung(1)=j*L;
						//Relativkoordinaten bestimmen
						VectorXd rel=Yn.col(x).head(row/2)-(Yn.col(xstrich).head(row/2)+Verschiebung);
						//Beschleunigung berrechnen
						A.col(x)+=f(rel, rc);

						if (x<xstrich)
						{
							//Potentielle Energie berechnen
							Epot += pot(rel, rc);
						}

					}
				}

			}
		}
	}
	return A;
}

//Funktion zur Berechnung von g
//Übergabewerte sind die Matrix Yn mit den Orts- und Geschwindigkeitsvektoren, die Größe L des Kastens,
//und der Vektor g in dem die Werte gespeichert werden
void Paarkorrelation(MatrixXd &Yn, double L, VectorXd &g) {
	int row= Yn.rows();
	int col= Yn.cols();
	int N=g.size()-1;
	// Binlänge bestimmen
	double dr=L/(2*N);
	//Teilchen Verschieben um virtuelle Teilchen zu generieren
	for (int xstrich = 0; xstrich < col; ++xstrich)
	{
		for (int x = 0; x < xstrich; ++x)
		{
			VectorXd Verschiebung=VectorXd::Zero(row/2);
			for (int i = -1; i <= 1; ++i)
			{
				Verschiebung(0)=i*L;
				for (int j = -1; j <= 1; ++j)
				{
					Verschiebung(1)=j*L;
					//Relativkoordinaten berechnen
					VectorXd rel=Yn.col(x).head(row/2)-(Yn.col(xstrich).head(row/2)+Verschiebung);

					double r =rel.norm();
					if (r<L/2)
					{
						//Anzahl der Paare zwischen r und r+dr bestimmen
						int l = (int) (r/dr);

						//g mit entsprechender Normierung, je Paar wird aufaddiert
						g(l)+=L*L/(col*col*M_PI*(-(l*dr)*(l*dr)+((l+1)*dr)*((l+1)*dr)));
					}
				}
			}
		}
	}
}

//Funktion zur Erstellung einer Ausgangsmatrix Y0 mit 16 Teilchen auf einem Gitter in einem Kasten der Länge L
// und zufälligen Geschwindigkeiten, wobei der Schwerpunktsimpuls 0 ist
MatrixXd Initialsierung(double L) {
	MatrixXd Y0(4,16);
	int j=0;
	//Teilchen auf Gitter erstellen
	for (int n = 0; n <= 3; ++n)
	{
		for (int m = 0; m <= 3; ++m)
		{
			Y0.col(j)(0)=L*(1+2.0*n)/8.0;
			Y0.col(j)(1)=L*(1+2.0*m)/8.0;
			Y0.col(j).tail(2)=VectorXd::Random(2);
			j+=1;
		}
	}
	//Schwerpunktsgeschwindigkeit berechnen
	VectorXd v_sum=VectorXd::Zero(2);
	for (int i = 0; i < 16; ++i)
	{
		v_sum +=Y0.col(i).tail(2);
	}
	//normierte Schwerpunktsgeschwindigkeit subtrahieren
	for (int i = 0; i < 16; ++i)
	{
		Y0.col(i).tail(2)-=v_sum/16;
	}

	//Startkonfiguration speichern
	ofstream bfile ("Data/start.txt", std::ofstream::out);
	for (int i = 0; i < 16; ++i)
	{
		bfile << Y0.col(i)(0) << "\t" << Y0.col(i)(1) << "\n";
	}
	return Y0;

}

//Funktion zur Ermittlung der Geschwindigkeit für Yn
//Übergeben werden die Matrix Y_{n+1}, die Matrix Y_n und Y_{n-1} sowie die Schrittweite h
void Velocity (MatrixXd &Yn_plus_1, MatrixXd &Yn, MatrixXd &Yn_minus_1, double &h){
	int row= Yn.rows();
	int col= Yn.cols();
	for (int i = 0; i < col; ++i)
	{
		Yn.col(i).tail(row/2)=(Yn_plus_1.col(i).head(row/2)-Yn_minus_1.col(i).head(row/2))/(2.0*h);
	}
}

//Funktion für die periodischen Randwerte
//Übergeben werden die Matrizen Y_n, Y_(n-1) und Y_(n-2) (3 Stück, um Sprünge in der Geschwindigkeit zu vermeiden)
// Sowie die Länge L der Box
void Periodic(MatrixXd &Yn,MatrixXd &Yn1,MatrixXd &Yn2, double L) {
	int row= Yn.rows();
	int col= Yn.cols();
	for (int j = 0; j < row/2; ++j)
	{
		for (int l = 0; l < col; ++l)
		{
			if(Yn.col(l)(j)<0 ||Yn.col(l)(j)>L)
			{
				Yn1.col(l)(j)=Yn1.col(l)(j)-L*floor(Yn.col(l)(j)/L);
				Yn2.col(l)(j)=Yn2.col(l)(j)-L*floor(Yn.col(l)(j)/L);
				Yn.col(l)(j)=Yn.col(l)(j)-L*floor(Yn.col(l)(j)/L);
			}
		}
	}
}

//Funktion zur Durchführung des Verlet Algorithmus
//Übergeben werden: die Funktion f zur Bestimmung der Kraft, die Matrizen Y_n und Y_(n-1), die Länge L des Kastens,
// die Schrittweite h, die Funktion pot zur Berechnung der potentiellen Energie und der double Epot zum Speichern der potentiellen Energie
MatrixXd Verlet (VectorXd (*f)(VectorXd &, double &),MatrixXd &Yn, MatrixXd &Yn_minus_1, double L, double h, double (*pot)(VectorXd &, double &) ,double &Epot) {
	int row= Yn.rows();
	int col= Yn.cols();
	MatrixXd Yn_plus_1= MatrixXd::Zero(row, col);

	MatrixXd a=accelaration(f,Yn, L, pot, Epot);

	for (int i = 0; i < col; ++i)
	{
		Yn_plus_1.col(i).head(row/2)=2*Yn.col(i).head(row/2)-Yn_minus_1.col(i).head(row/2)+a.col(i)*h*h;
	}
	return Yn_plus_1;
}

//Funktion zum skalieren der Temperatur
//Übergeben werden die Matrix Y mit den Orts und Geschwindigkeitsvektoren sowie die gewünschte Temperatur T0
void set_T(MatrixXd &Y, double &T0) {
	double T_momentan=0;
	int row= Y.rows();
	int col= Y.cols();
	//Momentane Geschwindigkeit berrechnen
	for (int i = 0; i < col; ++i)
	{
		T_momentan+=Y.col(i).tail(row/2).squaredNorm();
	}
	//Auf Teilchenzahl normieren
	T_momentan=T_momentan/((2.0*col-2.0));
	//Geschwindigkeit reskalieren, sodass T=T0
	for (int i = 0; i < col; ++i)
	{
		Y.col(i).tail(row/2)=Y.col(i).tail(row/2)*(sqrt(T0/T_momentan));
	}
}

//Funktion zur Äquilibrierung
//Übergeben werden ein string name für die Erstellung der txt Datei, die Funktion f für die Kraft, die Schrittweite h,
//Die Ausgangstemperatur T0, die Länge L des Kastens, die Anzahl der Schritte n, die Funktion pot zur Berechnung der potentiellen
//Energie, und ein String iskon (yes oder no) der festlegt, ob das isokinetische Thermostat verwendet wird
MatrixXd Aqu (string name,VectorXd (*f)(VectorXd &, double &), double h, double T0, double L, int n, double (*pot)(VectorXd &, double &), string isokin) {
	//Initialsierung aufrufen
	MatrixXd Y0=Initialsierung(L);

	//File zum Speichern einiger Zwischenkonfigurationen
	ofstream afile (give_name("Data/Werte_aqu_",name,".txt"), std::ofstream::out);

	int row= Y0.rows();
	int col= Y0.cols();

	//Vektor zum Speichern der potentiellen Energie erstellen
	VectorXd Epot=VectorXd::Zero(n);

	//Auf gewünschte Anfangstemperatur reskalieren
	set_T(Y0,T0);
	//Y_{-1} für Verlet berechnen
	MatrixXd Y_minus_1 =MatrixXd::Zero(row,col);
	MatrixXd a0=accelaration(f,Y0, L, pot, Epot(0));
	for (int i = 0; i < col; ++i)
	{
		Y_minus_1.col(i).head(row/2)=Y0.col(i).head(row/2)-Y0.col(i).tail(row/2)*h+0.5*a0.col(i)*h*h;
	}

	MatrixXd Yn=Y0;
	MatrixXd Yn_minus_1=Y_minus_1;

	//Matrix für den Schwerpunkt erstellen
	MatrixXd SP=MatrixXd::Zero(row/2,n+2);
	//Vektor zum Speichern der Temperatur erstellen
	VectorXd T=VectorXd::Zero(n);

	//Eigentliche Simulation starten
	for (int i = 0; i < n; ++i)
	{
		//Verlet Algorithmus aufrufen
		MatrixXd Yn_plus_1=Verlet(f,Yn, Yn_minus_1,L, h, pot, Epot(i));
		//Periodische Randbedingungen anwenden
		Periodic(Yn_plus_1,Yn, Yn_minus_1, L);
		//Geschwindigkeit ermitteln
		Velocity(Yn_plus_1, Yn, Yn_minus_1,h );
		//Prüfen, ob Thermostat verwendet werden soll und ggf, reskalieren
		if (isokin.compare("yes") == 0)
		{
			set_T(Yn, T0);
		}
		//Schwerpunktsgeschwindigkeit und Temperatur berechnen
		for (int j = 0; j < col; ++j)
		{
			SP.col(i)+=Yn.col(j).tail(row/2)/col;
			T(i)+= Yn.col(j).tail(row/2).squaredNorm()/(2.0*col-2.0);
		}
		//Jede 10ste Konfiguration speichern zur Visualisierung
		if (i % 10 == 0)
		{
			for (int j = 0; j < col; ++j)
				{
					afile << Yn.col(j)(0) << "\t" << Yn.col(j)(1) << "\n";
				}
		}
		//Matrizen umbennen
		Yn_minus_1=Yn;
		Yn=Yn_plus_1;
	}
	//Schwerpunktsgeschwindigkeit, Temperatur und potentielle Energie speichern
	ofstream dfile (give_name("Data/aqu_",name,".txt"), std::ofstream::out);
	for (int i = 0; i < n; ++i)
	{
		dfile << i*h << "\t" << SP.col(i)(0) << "\t" <<  SP.col(i)(1) << "\t" << T(i) << "\t" << Epot(i) << "\n";
	}
	return Yn_minus_1;
}

//Funktion zur eigentlichen Molekulardynamik Simulation
//Übergeben werden ein string für den Namen der txt datei, die Funktion f für die Kraft, die Schrittweite h,
//die Ausgangstemperatur T0, die Länge L des Kastens, die Schrittweite n, die Funktion pot zur Bestimmung der potentiellen Energie
//und ein string isokin (yes oder no) zur Festlegung, ob das isokinetische Thermostat verwendet wird
void MD (string name,VectorXd (*f)(VectorXd &, double &), double h, double T0, double L, int n,
double (*pot)(VectorXd &, double &), string isokin) {
	//Äquilibrierung durchführen
	MatrixXd Y0=Aqu(name,f,h,T0,L,1000, pot, isokin);
	int row= Y0.rows();
	int col= Y0.cols();

	//File zum Speichern einiger Zwischenkonfigurationen
	ofstream dfile (give_name("Data/Werte_",name,".txt"), std::ofstream::out);

	//Vektor zum Speichern der potentiellen Enerhie erstellen
	VectorXd Epot=VectorXd::Zero(n);
	//Y_{n-1} für Verlet berechnen
	MatrixXd Y_minus_1 =MatrixXd::Zero(row,col);
	MatrixXd a0=accelaration(f,Y0, L, pot, Epot(0));
	for (int i = 0; i < col; ++i)
	{
		Y_minus_1.col(i).head(row/2)=Y0.col(i).head(row/2)-Y0.col(i).tail(row/2)*h+0.5*a0.col(i)*h*h;
	}
	MatrixXd Yn=Y0;
	MatrixXd Yn_minus_1=Y_minus_1;
	//Vektor zum Speichern der Temperatur
	VectorXd T=VectorXd::Zero(n);
	//Vektor zum Speichern der Paarkorrelation
	VectorXd g=VectorXd::Zero(100);
	//Start der eigentlichen Simulation
	for (int i = 0; i < n; ++i)
	{
		//Aufruf des Verlet Algorithmus
		MatrixXd Yn_plus_1=Verlet(f,Yn, Yn_minus_1,L, h, pot, Epot(i));
		//Anwenden der periodischen Randbedingungen
		Periodic(Yn_plus_1,Yn, Yn_minus_1, L);
		//Bestimmung der Geschwindigkeiten
		Velocity(Yn_plus_1, Yn, Yn_minus_1,h );
		//Prüfe ob Thermostat verwendet werden soll und ggf. reskalieren
		if (isokin.compare("yes") == 0)
		{
			set_T(Yn, T0);
		}
		//Temperatur bestimmen
		for (int j = 0; j < col; ++j)
		{
			T(i)+= Yn.col(j).tail(row/2).squaredNorm()/(2.0*col-2.0);
		}
		//Paarkorrelation bestimmen
		Paarkorrelation(Yn, L, g);

		//Jede 10te Konfiguration speichern zur Visualisierung
		if (i % 10 == 0)
		{
			for (int j = 0; j < col; ++j)
				{
					dfile << Yn.col(j)(0) << "\t" << Yn.col(j)(1) << "\n";
				}
		}

		//Umbennenung für nächsten Schritt
		Yn_minus_1=Yn;
		Yn=Yn_plus_1;
	}

	//Temperatur mitteln
	double T_m=0;
	for (int i = 0; i < n; ++i)
	{
		T_m+= T(i)/n;

	}

	//Länge von g, also Anzahl der Bins bestimmen
	int gsize=g.size();
	//Paarkorrelation speichern
	ofstream afile (give_name("Data/g_",name,".txt"), std::ofstream::out);
	for (int i = 0; i < gsize; ++i)
	{
		afile << L/(2*gsize)*i << "\t" << g(i)/n << "\n";
	}
	ofstream bfile (give_name("Data/T_m_",name,".txt"), std::ofstream::out);
	bfile << T_m;
	//Temperatur und potentielle Energie speichern
	ofstream cfile (give_name("Data/T_",name,".txt"), std::ofstream::out);
	for (int i = 0; i < n; ++i)
	{
		cfile << i*h << "\t" << T(i) << "\t" << Epot(i) << "\n";

	}


}


int main()
{
	//Werte festsetzen
	double h=0.01; //Schrittweite
	double T1=1; //temperaturen
	double T2=0.01;
	double T3=10;

	double L=8.0; //Kastenlänge
	int n=1e5; //Anzahl der Schritte

	//Durchführung der Simulationen
	MD("1",Kraft_LJ,h,T1,L,n, Pot_LJ, "no");
	MD("0_01",Kraft_LJ,h,T2,L,n, Pot_LJ, "no");
	MD("10",Kraft_LJ,h,T3,L,n, Pot_LJ, "no");
	MD("isokin",Kraft_LJ,h,T2,L,n, Pot_LJ, "yes");



}
