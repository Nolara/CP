#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include <random>

using namespace std;
using namespace Eigen;

string give_name( const string& basename, string index, const string& ext) {
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

MatrixXi Initialise(int n, string random, mt19937 &generator) {

    MatrixXi Y=MatrixXi::Zero(n,n);
    uniform_int_distribution<int> start_distribution(0, 1); // liefert den Wert 0 oder 1
    if (random.compare("random") == 0) {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                Y(i,j)=start_distribution(generator)*2 - 1;
            }
        }
    } else {
        int start_direction =start_distribution(generator)*2 - 1;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                Y(i,j)=start_direction;
            }
        }
    }
    return Y;
}

void Aquilibrate(string name ,MatrixXi &Y,int n, int N, mt19937 &generator, uniform_real_distribution<double> &distribution , uniform_int_distribution<int> &pick_spin, double kbT) {
    double n_d=n;
    ofstream afile (give_name("Data/2_",name,".txt"), std::ofstream::out);
    VectorXd E=VectorXd::Zero(N+1);
	for (int k = 0; k < n; ++k)
	{
		for (int l = 0; l < n; ++l)
		{
			E(0)+=-Y(k,l)*(Y(k-1-n_d*floor((k-1)/n_d),l)+Y(k,l-1-n_d*floor((l-1)/n_d))+Y(k+1-n_d*floor((k+1)/n_d),l)+Y(k,l+1-n_d*floor((l+1)/n_d)));
		}
	}
	afile << Y << "\n";
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < n*n; ++j)
        {
            double random_number = distribution(generator);
            int x = pick_spin(generator);
            int y = pick_spin(generator);
            int delta_E=2*Y(x,y)*(Y(x-1-n_d*floor((x-1)/n_d),y)+Y(x,y-1-n_d*floor((y-1)/n_d))+Y(x+1-n_d*floor((x+1)/n_d),y)+Y(x,y+1-n_d*floor((y+1)/n_d)));

            //Überprüfen, ob der MC Move akzepiert wird
            if (delta_E<0) {
                Y(x,y)=-Y(x,y);
            } else if (random_number<exp(-delta_E/kbT)) {
                Y(x,y)=-Y(x,y);
            }

        }
        // if (N % 10 == 0)
        // {
        //     afile << Y << "\n";
        // }
        afile << Y << "\n";
        for (int k = 0; k < n; ++k)
        {
            for (int l = 0; l < n; ++l)
            {
                E(i+1)+=-Y(k,l)*(Y(k-1-n_d*floor((k-1)/n_d),l)+Y(k,l-1-n_d*floor((l-1)/n_d))+Y(k+1-n_d*floor((k+1)/n_d),l)+Y(k,l+1-n_d*floor((l+1)/n_d)));
            }
        }

    }
    ofstream bfile (give_name("Data/2_E_",name,".txt"), std::ofstream::out);
    bfile << E/(n*n);
    bfile.close();
    afile.close();
}

void measure (string name ,MatrixXi &Y,int n, int N, mt19937 &generator, uniform_real_distribution<double> &distribution , uniform_int_distribution<int> &pick_spin, double kbT) {
    double n_d=n;
    VectorXd m=VectorXd::Zero(N+1);
    VectorXd E=VectorXd::Zero(N+1);
	for (int k = 0; k < n; ++k)
	{
		for (int l = 0; l < n; ++l)
		{
			m(0)+=Y(k,l);
			E(0)+=-Y(k,l)*(Y(k-1-n_d*floor((k-1)/n_d),l)+Y(k,l-1-n_d*floor((l-1)/n_d))+Y(k+1-n_d*floor((k+1)/n_d),l)+Y(k,l+1-n_d*floor((l+1)/n_d)));
		}
	}
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < n*n; ++j)
        {
            double random_number = distribution(generator);
            int x = pick_spin(generator);
            int y = pick_spin(generator);
            int delta_E=2*Y(x,y)*(Y(x-1-n_d*floor((x-1)/n_d),y)+Y(x,y-1-n_d*floor((y-1)/n_d))+Y(x+1-n_d*floor((x+1)/n_d),y)+Y(x,y+1-n_d*floor((y+1)/n_d)));

            //Überprüfen, ob der MC Move akzepiert wird
            if (delta_E<0) {
                Y(x,y)=-Y(x,y);
            } else if (random_number<exp(-delta_E/kbT)) {
                Y(x,y)=-Y(x,y);
            }

        }
        for (int k = 0; k < n; ++k)
        {
            for (int l = 0; l < n; ++l)
            {
                m(i+1)+=Y(k,l);
                E(i+1)+=-Y(k,l)*(Y(k-1-n_d*floor((k-1)/n_d),l)+Y(k,l-1-n_d*floor((l-1)/n_d))+Y(k+1-n_d*floor((k+1)/n_d),l)+Y(k,l+1-n_d*floor((l+1)/n_d)));
            }
        }

    }
	double H_expectation=0;
	double H_squared=0;
	for (int k = 0; k < N+1; ++k)
	{
		H_expectation+=E(k);
		H_squared+=E(k)*E(k);
	}
	H_expectation=H_expectation/N;
	H_squared=H_squared/N;
	double c=0;
	c=(H_squared-H_expectation*H_expectation)/(kbT*kbT*n*n);
	m=m/(n*n);
    ofstream bfile (give_name("Data/2_m_",name,".txt"), std::ofstream::out);
    bfile << "# C=" << c << "\n";
    bfile << m;
	bfile.close();
}


int main(){
    std::random_device rd;
	std::mt19937 generator(rd());

    int n=100;
    double kbT1=1;
    double kbT1_5=1.5;
    double kbT2_27=2.27;
    double kbT3=3;


    uniform_int_distribution<int> pick_spin(0, n-1);
    uniform_real_distribution<double> distribution(0,1);

    MatrixXi Y_random3=Initialise(n, "random", generator);
    MatrixXi Y_ordered1=Initialise(n, "ordered", generator);
    MatrixXi Y_random1=Initialise(n, "random", generator);
    MatrixXi Y_ordered3=Initialise(n, "ordered", generator);

    Aquilibrate("3_random",Y_random3,n,100,generator,distribution, pick_spin, kbT3);
    Aquilibrate("1_random",Y_random1,n,100,generator,distribution, pick_spin, kbT1);
    Aquilibrate("1",Y_ordered1,n,100,generator,distribution, pick_spin, kbT1);
    Aquilibrate("3",Y_ordered3,n,100,generator,distribution, pick_spin, kbT3);

    MatrixXi Y_random1_5=Initialise(n, "random", generator);
    MatrixXi Y_ordered1_5=Initialise(n, "ordered", generator);
    MatrixXi Y_random2_27=Initialise(n, "random", generator);
    MatrixXi Y_ordered2_27=Initialise(n, "ordered", generator);

    Aquilibrate("1_5_random",Y_random1_5,n,100,generator,distribution, pick_spin, kbT1_5);
    Aquilibrate("2_27_random",Y_random2_27,n,100,generator,distribution, pick_spin, kbT2_27);
    Aquilibrate("1_5",Y_ordered1_5,n,100,generator,distribution, pick_spin, kbT1_5);
    Aquilibrate("2_27",Y_ordered2_27,n,100,generator,distribution, pick_spin, kbT2_27);

	measure("1_5",Y_ordered1_5,n,1000,generator,distribution, pick_spin, kbT1_5);
	measure("2_27",Y_ordered2_27,n,1000,generator,distribution, pick_spin, kbT2_27);
	measure("3",Y_ordered3,n,1000,generator,distribution, pick_spin, kbT3);

    measure("1_5_random",Y_random1_5,n,1000,generator,distribution, pick_spin, kbT1_5);
    measure("2_27_random",Y_random2_27,n,1000,generator,distribution, pick_spin, kbT2_27);
    measure("3_random",Y_random3,n,1000,generator,distribution, pick_spin, kbT3);


}
