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
	int N(3), k(1), k_max(10000);
	double eps(0.000000001), a(0.1);
	MatrixXd C(N,N);
	SparseMatrix<double> Id(N,N), A(N,N), B(N,N), M(N,N), N_J(N,N), V(N,N), H(N,N);
	SparseVector<double> b(N), sol0(N), sol(N), r(N);
	int userchoicemethode;
	string results;


//NE PAS SUPPRIMER LES COMMENTAIRES D'INITIALISATION SVP	 

	// for (int i=0 ; i<A.rows() ; ++i)
	// {
	//   for (int j=0 ; j<A.cols() ; ++j)
	//   {
	// 	A.coeffRef(i,j)=1;
	//   }
	// }

//	Définition des matrices à utiliser globalement
	Id.setIdentity();                          //Matrice Identité
	C=MatrixXd::Random(N,N);     //Matrice random B
	B = C.sparseView();
	A=a*Id+B.transpose()*B;
	//Matrice A

  // MatrixXd Unit(N,N);
	// // C est une matrice de réels aléatoires entre -1 et 1
	// C=MatrixXd::Random(N,N);
	// // Unit est la matrice dont les coefficients sont tous égaux à 1
	// Unit=MatrixXd::Constant(N,N,1.);
	// // La matrice C est désormais une matrice de réels aléatoires entre 0 et 1
	// C=1./2.*(Unit+C);
	// // B est la représentation creuse de C
	// B=C.sparseView();
	// // définition de la matrice identité
//	Id.setIdentity();
	// Initialisation du vecteur b
	//A= 3*N*Id+B;


	for (int i=0 ; i<sol0.rows() ; i++)
	{
		sol0.coeffRef(i)=1;
		b.coeffRef(i)=1;                 //Définir un valeur de sol0
	}

	r=b-A*sol0;                                //Initialisation de r


	cout << "------------------------------------" << endl;
	cout << "Choississez la méthode de résolution : " << endl;
	cout << "1) Jacobi"<< endl;
	cout << "2) Gradient Pas Optimal" << endl;
	cout << "3) Résidu Minimum" << endl;
	//cout << "4) GMRes" << endl;
	cin >> userchoicemethode;

	MethodeRes* methode(0);

	switch(userchoicemethode)
	{

		case 1: //Jacobi
		methode = new Jacobi( M, N_J);
		results = "solution_Jacobi.txt";    // Nom du fichier solution
		break;


		case 2: //GPO
		methode = new GPO();
		results = "solution_GPO.txt";           // Nom du fichier solution
		break;


		case 3: //Résidu
		methode = new Residu();
		results = "solution_Residu.txt";    // Nom du fichier solution
		break;

		case 4: //GMRes
		methode = new GMRes(V, H);
		results = "solution_GMRes.txt";    // Nom du fichier solution
		break;

		default:
		cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
		exit(0);
	}

	// Initialisations
	methode->Initialisation(b,A,sol0,r,results,methode);

	// On sauvegarde la solution
	methode->SaveSolution(0, r);

	while (r.norm()>eps && k<=k_max)
	{
		methode->calcul_sol(r);   //Appel de la fonction solution
		methode->SaveSolution(k,r);

		cout << "k=" << k << "\n"<< endl;
		cout << "=======================" << endl;
		//cout << "rbouc" <<r.norm() << endl;
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
