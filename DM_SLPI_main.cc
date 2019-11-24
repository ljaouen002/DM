#include "Fonctions.h"
#include "MethodeRes.h"
#include <string>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Eigen;

int main()
{
	int N(10);
  	SparseMatrix Id(N,N), A(N,N), D(N,N), E(N,N), F(N,N), B(N,N);
	SparseVector sol(N);
  	double alpha(0.1);

  	A.setZero(A.rows(),A.cols()); E.setZero(A.rows(),A.cols()); F.setZero(A.rows(),A.cols());
  
  	Id = MatrixXd::Identity(10,10);    //Matrice Identité
	B=MatrixXd::Random(10,10);        //Matrice random B
  	A=alpha*I+B.transpose()*B;        //Matrice A
  
  	D=diag(A);   // Diagonale de A
	for (int i=0 ; i<A.rows() ; i++)
	{
		for (int j=0 ; j<A.cols() ; j++)
		{
			E=-A(i,j+1);     // Partie triangulaire supérieure de A
  			F=-A(j+1,i);     // Partie triangulaire inférieure de A
		}
	}
	
	return 0;
}
