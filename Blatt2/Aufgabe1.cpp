#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>


using namespace std;
using namespace Eigen;


int loadData(Eigen::MatrixXd &mat, std::string filename, size_t row, size_t col)
{

	std::ifstream input;
	input.open(filename, std::fstream::in | std::fstream::binary);
	if (input.is_open())
	{
		std::vector<unsigned char> vec((std::istreambuf_iterator<char>(input)),
				std::istreambuf_iterator<char>() );
		mat= Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic,Eigen::Dynamic>>(vec.data(),row,col).template cast<double>();
		input.close();
		return 1;
	}
	return 0;

}
int storeData(Eigen::MatrixXd &mat, std::string filename)
{
	Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic> pic= mat.template cast<char>();
	std::ofstream output;
	output.open(filename,std::fstream::trunc);
	if (output.is_open())
	{
		output.write(pic.data(),pic.size());
		output.close();
		return 1;
	}
	return 0;
}

string name( const string& basename, float index, const string& ext)
{
	ostringstream name;
	name << basename << index << ext;
	return name.str();
}

void Reduktion(int k, MatrixXd Singular, MatrixXd U, MatrixXd V, int n)
{
    MatrixXd Reduktion(n, n);

    for (float i = 0; i < k; i++)
    {
    	Reduktion += Singular(i) * U.col(i) * V.row(i);
    }

    ofstream file;
    file.open (name("Data/Reduktion_", k, ".txt"), ofstream::out);
    file << Reduktion.transpose();
    file.close();
    Reduktion.resize(0,0);
}

int main()
{
    int size = 512;
    MatrixXd M(size, size);
    loadData(M, "Bild", size, size);
    MatrixXd Original(size, size);
    Original = M.transpose();

    BDCSVD<MatrixXd> svd(M, ComputeFullV | ComputeFullU );

    MatrixXd Singular(size, size);
    Singular=svd.singularValues();

    MatrixXd U(size, size);
    U=svd.matrixU();

    MatrixXd V(size, size);
    V=svd.matrixV().transpose();


    ofstream afile ("Data/one_original.txt", std::ofstream::out); // Erstelle txt Datei für Aufgabenteil a
    afile << "# Aufgabenteil b:" << "\n" <<  "\n";
    afile << Original << endl;
    afile.close();

    Reduktion(10, Singular,U, V, size);
    Reduktion(20, Singular,U, V, size);
    Reduktion(50, Singular,U, V, size);

}
