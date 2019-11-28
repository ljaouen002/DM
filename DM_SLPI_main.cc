//#include "Fonctions.h"
#include "MethodeRes.h"
#include <string>
#include <iostream>
#include <cmath>
#include "Dense"
#include "Sparse"

using namespace std;
using namespace Eigen;

int main()
{
	int nb_iterations(100);
	int N(10), k(0), k_max(100);
	double eps(0.0001);
	SparseMatrix<double> Id(N,N), A(N,N), D(N,N), E(N,N), F(N,N), B(N,N);
	SparseVector<double> sol(N), sol0(N), b(N);
	double alpha_sujet(0.1);
	SparseVector<double> r;
	int userchoicemethode;
	string results;

	// Définition des matrices à utiliser globalement

	Id=MatrixXd::Identity(N,N);    //Matrice Identité
	B.setRandom(B.rows(),B.cols());        //Matrice random B
	//Lisa: Ce que je comprends du sujet, c'est que la matrice  ne doit contenir que des valeurs égale à 0 ou à 1, pas entre les deux.
	A=alpha_sujet*Id+B.transpose()*B;        //Matrice A


	cout << "------------------------------------" << endl;
	cout << "Choississez la méthode de résolution : " << endl;
	cout << "1) Jacobi"<< endl;
	cout << "2) Gradient Pas Optimal" << endl;
	cout << "3) Résidu Minimum" << endl;
	cout << "4) GMRes" << endl;
	cin >> userchoicemethode;
	//userchoicemethode=1;

	MethodeRes* methode(0);

	switch(userchoicemethode)
	{

		case 1: //Jacobi
			methode = new Jacobi();
			// Nom du fichier solution
			results = "solution_Jacobi.txt";
		break;


		case 2: //GPO
			double alpha; //Attention, alpha ici est le coefficient de descente, différent du alpha précédent de l'énoncé
			methode = new GPO(alpha);
			results = "solution_GPO";
		break;


		case 3: //Résidu
			double alpha;
			methode = new Residu(alpha);
			results = "solution_Residu.txt";
		break;

		case 4: //GMRes
			//double beta;
			//methode = new GMRes(beta);
			//results = "solution_GMRes.txt";
		break;

		default:
		cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
		exit(0);
	}

	// Initialisations
	methode->Initialisation(b,A,sol0,methode);

	// On sauvegarde la solution
	methode->SaveSolution(0);

	//Faire une boucle pour trace la norme de _r en fonction de nb_iterations ????
	while (r.norm()>eps || k<=k_max)
	{
		methode->calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);   //Appel de la fonction solution
		methode->SaveSolution(k);
		k+=1;
	}
	//cout << _sol << endl;
	if (k>k_max)
	{
		cout << "Tolérance non atteinte :" << r.norm() << endl;
	}


	delete methode;
	return 0;
}
