#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using namespace Eigen;

string give_name( const string& basename, string index, const string& ext) {
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

MatrixXd Initialise(double Delta) {
    int N=18/Delta;
    MatrixXd r0=MatrixXd::Zero(N,2);
    for (int i = 0; i < 1/Delta; ++i)
    {
        r0.row(i) << i*Delta, 0;
    }
    for (int i = 1/Delta; i < 3/Delta; ++i)
    {
        r0.row(i) << 1, i*Delta-1;
    }
    for (int i = 3/Delta; i < 5/Delta; ++i)
    {
        r0.row(i) << i*Delta-2, 2;
    }
    for (int i = 5/Delta; i < 7/Delta; ++i)
    {
        r0.row(i) << 3, 2-i*Delta+5;
    }
    for (int i = 7/Delta; i < 8/Delta; ++i)
    {
        r0.row(i) << i*Delta-4, 0;
    }
    for (int i = 8/Delta; i < 11/Delta; ++i)
    {
        r0.row(i) << 4 , i*Delta-8;
    }
    for (int i = 11/Delta; i < 15/Delta; ++i)
    {
        r0.row(i) << 4-i*Delta+11 , 3;
    }
    for (int i = 15/Delta; i < 18/Delta; ++i)
    {
        r0.row(i) << 0 , 3-i*Delta+15;
    }
    return r0;

}

void Simulated_Annealing (MatrixXd &r_i, VectorXi &permutation, double T_start, double T_end, double d, int s, string name) {
    double T=T_start;
    int N=r_i.rows();
    double N_d=N;
    MatrixXd r_i_plus_1=r_i;

    std::random_device rd;
	std::mt19937 generator(rd());

    int n = (int) (log(T_end/T_start)/log(d));
    VectorXd L=VectorXd::Zero(n+2);


    for (int i = 0; i < N-1; ++i)
    {
        L(0)+=(r_i.row(i)-r_i.row(i+1)).norm();
    }
    L(0)+=(r_i.row(N-1)-r_i.row(0)).norm();

    uniform_int_distribution<int> pick_space(0,N-1);
    uniform_real_distribution<double> distribution(0,1);
    int index=1;
    while (T>T_end)
    {

        for (int i = 0; i < s; ++i)
        {
            int pick_space_1 = pick_space(generator);
            int pick_space_2 = pick_space(generator);
            if (pick_space_1!=pick_space_2)
            {
                double L_now=(r_i.row(pick_space_1)-r_i.row(pick_space_2-1-N_d*floor((pick_space_2-1)/N_d))).norm()+(r_i.row(pick_space_1)-r_i.row(pick_space_2+1-N_d*floor((pick_space_2+1)/N_d))).norm()+(r_i.row(pick_space_2)-r_i.row(pick_space_1-1-N_d*floor((pick_space_1-1)/N_d))).norm()+(r_i.row(pick_space_2)-r_i.row(pick_space_1+1-N_d*floor((pick_space_1+1)/N_d))).norm();
                double L_previous=(r_i.row(pick_space_1)-r_i.row(pick_space_1-1-N_d*floor((pick_space_1-1)/N_d))).norm()+(r_i.row(pick_space_1)-r_i.row(pick_space_1+1-N_d*floor((pick_space_1+1)/N_d))).norm()+(r_i.row(pick_space_2)-r_i.row(pick_space_2-1-N_d*floor((pick_space_2-1)/N_d))).norm()+(r_i.row(pick_space_2)-r_i.row(pick_space_2+1-N_d*floor((pick_space_2+1)/N_d))).norm();
                double delta_L=L_now-L_previous;
                double random_number=distribution(generator);
                if (delta_L<0 || exp(-delta_L/T)>random_number)
                {
                    int temp=permutation(pick_space_1);
                    VectorXd temp_vec=r_i.row(pick_space_1);
                    permutation(pick_space_1)=permutation(pick_space_2);
                    permutation(pick_space_2)=temp;
                    r_i.row(pick_space_1)=r_i.row(pick_space_2);
                    r_i.row(pick_space_2)=temp_vec;

                }
            }

        }

        for (int j = 0; j < N-1; ++j)
        {
            L(index)+=(r_i.row(j)-r_i.row(j+1)).norm();
        }
        L(index)+=(r_i.row(N-1)-r_i.row(0)).norm();
        T=T*d;
        index+=1;
    }
    ofstream cfile (give_name("Data/3_",name,".txt"), std::ofstream::out);
    cfile << L;
    ofstream dfile (give_name("Data/3_final_",name,".txt"), std::ofstream::out);
    dfile << r_i;
    ofstream efile (give_name("Data/3_permutation_",name,".txt"), std::ofstream::out);
    efile << permutation;
}

int main()
{
    MatrixXd r0=Initialise(0.2);
    ofstream afile ("Data/3_start.txt", std::ofstream::out);
    afile << r0;

    std::vector<int> permutation_std;
    for (int i = 0; i < r0.rows(); ++i)
    {
        permutation_std.push_back(i);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle (permutation_std.begin(), permutation_std.end(), std::default_random_engine(seed));
    Eigen::VectorXi permutation = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(permutation_std.data(), permutation_std.size());

    MatrixXd r0_random=MatrixXd::Zero(r0.rows(),2);
    for (int i = 0; i < r0.rows(); ++i)
    {
        r0_random.row(i)=r0.row(permutation(i));
    }
    ofstream bfile ("Data/3_start_random.txt", std::ofstream::out);
    bfile << r0_random;

    VectorXd d=VectorXd::Zero(3);
    d << 0.9,0.99,0.999;
    VectorXi s=VectorXi::Zero(5);
    s << 10,100,1e3,1e4,1e5;


    double T_start=10;
    double T_end=1e-2;
    int index=0;
    for (int i = 0; i <3; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            auto name = std::to_string(index);
            MatrixXd r=r0_random;
            VectorXi p=permutation;
            Simulated_Annealing(r, p, T_start, T_end, d(i), s(j), name);
            index+=1;
        }
    }





}
