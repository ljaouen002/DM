#include "MethodeRes.h"
#include <string>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "Dense"
#include "Sparse"

using namespace std;
using namespace Eigen;

int main()
{
	int N, k(1), k_max(300);
	double eps(0.001), a(0.1);
	MatrixXd C, Unit;
	SparseMatrix<double> Id, A, B, M, N_J, V, H, M_precon;
	SparseVector<double> b, sol0, sol, r, q;
	int userchoicemethode;
	int userchoicematrice;
	int precondi;
	int m_GMRes(10);
	string results, name_Matrix;

	// Choix de la matrice à utiliser
	cout << "------------------------------------" << endl;
	cout << "Choississez la matrice à utiliser : " << endl;
	cout << "1) Matrice avec seulement des 1"<< endl;
	cout << "2) Matrice du DM" << endl;
	cout << "3) Matrice marraine" << endl;
	cout << "4) Matrice BCSSTK18" << endl;
	cout << "5) Matrice FS_541_4" << endl;
	cin >> userchoicematrice;

	MethodeRes* matrice(0);

	switch(userchoicematrice)
	{
		case 1: //Matrice avec seulement des 1
		N=3;
		A.resize(N,N);

		for (int i=0 ; i<A.rows() ; ++i)
		{
			for (int j=0 ; j<A.cols() ; ++j)
			{
				A.coeffRef(i,j)=2;
			}
		}

		// Définition des vecteurs sol0 et b
		sol0.resize(N) ; b.resize(N) ; r.resize(N);
		for (int i=0 ; i<sol0.rows() ; i++)
		{
			sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
			b.coeffRef(i)=1.;          //Définir un valeur de b
		}

		r=b-A*sol0;                 //Initialisation de r
		break;


		case 2: //Matrice du DM
		N=100;
		Id.resize(N,N) ; C.resize(N,N) ; B.resize(N,N) ; A.resize(N,N);

		Id.setIdentity();              // Matrice Identité

		C = MatrixXd::Random(N,N);     // Matrice random C dense

		for (int i=0 ; i<A.rows() ; ++i)
		{
			for (int j=0 ; j<A.cols() ; ++j)
			{
				C(i,j)=abs(C(i,j));
			}
		}

		B = C.sparseView();            // Matrice random B sparse
		A = a*Id+B.transpose()*B;      // Matrice A
		//Note= pour que Jacobi converve, il faut multiplier a*n2*Id parce que sinon pas strcitement dominante


		// Définition des vecteurs sol0 et b
		sol0.resize(N) ; b.resize(N) ; r.resize(N);
		for (int i=0 ; i<sol0.rows() ; i++)
		{
			sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
			b.coeffRef(i)=1.;          //Définir un valeur de b
		}

		r=b-A*sol0;                 //Initialisation de r
		break;


		case 3:
		N=100;
		Id.resize(N,N) ; C.resize(N,N) ; B.resize(N,N) ; A.resize(N,N); Unit.resize(N,N);

		// C est une matrice de réels aléatoires entre -1 et 1
		C = MatrixXd::Random(N,N);
		// Unit est la matrice dont les coefficients sont tous égaux à 1
		Unit = MatrixXd::Constant(N,N,1.);
		// La matrice C est désormais une matrice de réels aléatoires entre 0 et 1
		C = 1./2.*(Unit+C);
		// B est la représentation creuse de C
		B = C.sparseView();
		// Définition de la matrice identité
		Id.setIdentity();
		// Initialisation du vecteur b
		A = 3*N*Id+B;

		// Définition des vecteurs sol0 et b
		sol0.resize(N) ; b.resize(N) ; r.resize(N);
		for (int i=0 ; i<sol0.rows() ; i++)
		{
			sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
			b.coeffRef(i)=1.;          //Définir un valeur de b
		}

		r=b-A*sol0;                 //Initialisation de r
		break;

		case 4: //BCSSTK18
		name_Matrix = "BCSSTK18";
		matrice->InitialisationMat(name_Matrix,A,b,sol0,r);
		break;


		case 5: //FS_541_4
		name_Matrix = "FS_541_4";
		matrice->InitialisationMat(name_Matrix,A,b,sol0,r);
		break;


		default:
		cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
		exit(0);
	}


	// Choix de la méthode de résolution
	cout << "------------------------------------" << endl;
	cout << "Choississez la méthode de résolution : " << endl;
	cout << "1) Jacobi"<< endl;
	cout << "2) Gradient Pas Optimal" << endl;
	cout << "3) Résidu Minimum" << endl;
	cout << "4) GMRes" << endl;
	cin >> userchoicemethode;

	MethodeRes* methode(0);

	switch(userchoicemethode)
	{

		case 1: //Jacobi
		M.resize(N,N);
		N_J.resize(N,N);
		methode = new Jacobi( M, N_J);
		results = "solution_Jacobi.txt";    // Nom du fichier solution
		break;


		case 2: //GPO
		methode = new GPO();
		results = "solution_GPO.txt";           // Nom du fichier solution
		break;


		case 3: //Résidu

		cout << "------------------------------------" << endl;
		cout << "Avec préconditionnement ? " << endl;
		cout << "0) Sans préconditionnement" << endl;
		cout << "1) à gauche par Jacobi"<< endl;
		cout << "2) à gauche par SGS" << endl;
		cout << "3) à droite par Jacobi"<< endl;
		cout << "4) à droite par SGS" << endl;
		cout << "5) auto-préconditon" << endl;
		cin >> precondi;

		if (precondi == 0)
		{
			methode = new Residu();
			results = "solution_Residu.txt";    // Nom du fichier solution
		}
		else if (precondi == 1)
		{
			methode = new Residu_Precondi_gauche(M_precon, q, precondi);
			results = "solution_Residu_Precondi_G_Jacobi.txt";    // Nom du fichier solution
		}
		else if (precondi == 2)
		{
			methode = new Residu_Precondi_gauche(M_precon, q, precondi);
			results = "solution_Residu_Precondi_G_SGS.txt";    // Nom du fichier solution
		}
		else if (precondi == 3)
		{
			methode = new Residu_Precondi_droite(M_precon, precondi);
			results = "solution_Residu_Precondi_D_Jacobi.txt";    // Nom du fichier solution
		}
		else if (precondi == 4)
		{
			methode = new Residu_Precondi_droite(M_precon, precondi);
			results = "solution_Residu_Precondi_D_SGS.txt";    // Nom du fichier solution
		}
		else if (precondi == 5)
		{
			methode = new Residu_Precondi_auto(M_precon);
			results = "solution_Residu_Precondi_Auto.txt";    // Nom du fichier solution
		}
		break;

		case 4: //GMRes
		V.resize(N,m_GMRes+1);
		H.resize(m_GMRes+1,m_GMRes);

		methode = new GMRes(V, H, m_GMRes);
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
		cout << "norme" <<r.norm() << endl;
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
