#include "Fonctions.h"
#include "MethodeRes.h"
#include <string>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Eigen;

int main()
{
	int nb_iterations(100);
	int N(10), k(0), k_max(100);
	double eps(0.0001);
  	SparseMatrix Id(N,N), A(N,N), D(N,N), E(N,N), F(N,N), B(N,N);
	SparseVector sol(N);
  	double alpha(0.1);
	
	// Définition d'un pointeur de MethodeRes
	MethodeRes* methode = new Jacobi(); 
	// Nom du fichier solution
	string results = "solution.txt";

	// Définition des matrices à utiliser
  	A.setZero(A.rows(),A.cols()); E.setZero(A.rows(),A.cols()); F.setZero(A.rows(),A.cols());
  
  	Id = MatrixXd::Identity(10,10);    //Matrice Identité
	B=MatrixXd::Random(10,10);        //Matrice random B
  	A=alpha*I+B.transpose()*B;        //Matrice A
  
  	D=diag(A);   // Diagonale de A
	for (int i=0 ; i<A.rows() ; i++)
	{
		for (int j=0 ; j<A.cols() ; j++)
		{
			if (i<j)
			{
				E=-A(i,j);     // Partie triangulaire supérieure de A
			}
			else if (i>j)
			{
				F=-A(i,j);     // Partie triangulaire inférieure de A
			}
		}
	}
	
	// Initialisation
	methode->Initialisation(b,A,sol0,methode);
	// On sauvegarde la solution
	methode->SaveSolution(nb_iterations);
	
	//Faire une boucle pour trace la norme de _r en fonction de nb_iterations ????
	
	// Algorithme de Jacobi
	 while (_r.norm()>eps || k<=k_max)
 	 {
   		k+=1
		methode->calcul_sol(Eigen::SparseVector b, Eigen::SparseMatrix A);   //Appel de la fonction Jacobi
	 }
	 //cout << _sol << endl;
  	if (k>k_max)
  	{
    	cout << "Tolérance non atteinte :" << _r.norm() << endl;
  	}
	
		
	delete methode;
	return 0;
}
