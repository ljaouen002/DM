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
	int N(10), k(0), k_max(100);
	double eps(0.0001), a(0.1);
	MatrixXd Id(N,N), A(N,N), B(N,N);
	MatrixXd D(N,N), E(N,N), F(N,N), M(N,N), N_J(N,N);
	VectorXd b(N), sol0(N), sol(N), r(N);
	int userchoicemethode;
	string results;

	// Définition des matrices à utiliser globalement
	Id=MatrixXd::Identity(N,N);                //Matrice Identité
	B=MatrixXd::Random(B.rows(),B.cols());     //Matrice random B
	A=a*Id+B.transpose()*B;                    //Matrice A
	sol0=VectorXd::Random(sol0.rows());        //Définir un valeur de sol0
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


			methode = new Jacobi(D, F, E, M, N_J);
			results = "solution_Jacobi.txt";    // Nom du fichier solution
		break;


		case 2: //GPO
			methode = new GPO();
			results = "solution_GPO";           // Nom du fichier solution
		break;


		case 3: //Résidu
			methode = new Residu();
			results = "solution_Residu.txt";    // Nom du fichier solution
		break;

		/*case 4: //GMRes
			methode = new GMRes(beta);
			results = "solution_GMRes.txt";    // Nom du fichier solution
		break;*/

		default:
			cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
			exit(0);
	}

	// Initialisations
	methode->Initialisation(b,A,sol0,r,results,methode);

	// On sauvegarde la solution
	methode->SaveSolution(0);

	while (r.norm()>eps && k<=k_max)
	{
		methode->calcul_sol();   //Appel de la fonction solution
		methode->SaveSolution(k);
		cout << "k=" << k << endl;
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
