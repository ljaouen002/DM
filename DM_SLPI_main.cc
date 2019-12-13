#include "MethodeRes.h"
#include "ReadMatrix.h"
#include <string>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "Dense"
#include "Sparse"
#include <chrono>

using namespace std;
using namespace Eigen;

int main()
{
	int N, k(1), k_max;
	double eps, a(0.1);
	MatrixXd C;
	SparseMatrix<double> Id, A, B, M, N_J, V, H, M_precon, D, D1, E, F;
	SparseVector<double> b, sol0, sol, r, q;
	int userchoicemethode;
	int userchoicematrice;
	int precondi;
	int m_GMRes(10);
	string results, name_Matrix;


	// Choix de la matrice à utiliser
	cout << "------------------------------------" << endl;
	cout << "Choississez la matrice à utiliser : " << endl;
	cout << "1) Matrice A random" << endl;
	cout << "2) Matrice BCSSTK18" << endl;
	cout << "3) Matrice FS_541_4" << endl;
	cin >> userchoicematrice;

	MethodeRes* matrice(0);

	switch(userchoicematrice)
	{
		case 1: //Première matrice du DM

		cout << "------------------------------------" << endl;
		cout << "Quelle taille de matrice voulez vous? " << endl;
		cin >> N;
		Id.resize(N,N) ; C.resize(N,N) ; B.resize(N,N) ; A.resize(N,N);

		Id.setIdentity();              				// Matrice Identité

		C = MatrixXd::Random(N,N);    			 	// Matrice random C dense

		for (int i=0 ; i<A.rows() ; ++i)
		{
			for (int j=0 ; j<A.cols() ; ++j)
			{
				C(i,j)=abs(C(i,j));								//C est composé de terme compris entre 0 et 1
			}
		}

		B = C.sparseView();            				// Matrice random B sparse
		A = a*Id+B.transpose()*B;     			 	// Matrice A
//		A = a*N*N*N*Id+B.transpose()*B;     // Matrice A pour l'étude de Jacobi

		// Définition des vecteurs sol0 et b
		sol0.resize(N) ; b.resize(N) ; r.resize(N);
		for (int i=0 ; i<sol0.rows() ; i++)
		{
			sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
			b.coeffRef(i)=1.;          //Définir un valeur de b
		}

		r=b-A*sol0;                 //Initialisation de r
		break;


		case 2: //BCSSTK18
			InitialisationMatrixA(N,"bcsstk18",A); 
			sol0.resize(N);
			b.resize(N);
			for (int i=0 ; i<sol0.rows() ; i++)
			{
				sol0.coeffRef(i)=1;        //Définir un valeur de sol0
				b.coeffRef(i)=1.;          //Définir un valeur de b
			}
		break;


		case 3: //FS_541_4
			InitialisationMatrixA(N,"fs_541_4",A);
			sol0.resize(N);
			b.resize(N);
			for (int i=0 ; i<sol0.rows() ; i++)
			{
				sol0.coeffRef(i)=1;        //Définir un valeur de sol0
				b.coeffRef(i)=1.;          //Définir un valeur de b
			}
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
		results = "solution_Jacobi_BC.txt";
		break;


		case 2: //GPO
		methode = new GPO();
		results = "solution_GPO_Compa.txt";
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
			results = "solution_Residu.txt";
		}
		else if (precondi == 1)
		{
			methode = new Residu_Precondi_gauche(M_precon, F, E, D, D1, q, precondi);
			results = "solution_Residu_Precondi_G_Jacobi.txt";
		}
		else if (precondi == 2)
		{
			methode = new Residu_Precondi_gauche(M_precon, F, E, D, D1, q, precondi);
			results = "solution_Residu_Precondi_G_SGS.txt";
		}
		else if (precondi == 3)
		{
			methode = new Residu_Precondi_droite(M_precon, F, E, D, D1, precondi);
			results = "solution_Residu_Precondi_D_Jacobi.txt";
		}
		else if (precondi == 4)
		{
			methode = new Residu_Precondi_droite(M_precon, F, E, D, D1, precondi);
			results = "solution_Residu_Precondi_D_SGS.txt";
		}
		else if (precondi == 5)
		{
			methode = new Residu_Precondi_auto();
			results = "solution_Residu_Precondi_Auto.txt";
		}
		break;

		case 4: //GMRes
		V.resize(N,m_GMRes+1);
		H.resize(m_GMRes+1,m_GMRes);

		methode = new GMRes(V, H, m_GMRes);
		results = "solution_GMRes_Compa.txt";
		break;

		default:
		cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
		exit(0);
	}

	//Choix du nombre d'itération maximum
	cout << "------------------------------------" << endl;
	cout << "Combien d'itération voulez vous? " << endl;
	cin >> k_max;

	//Choix du critère de convergence
	cout << "------------------------------------" << endl;
	cout << "Quel seuil de tolérance voulez vous? " << endl;
	cin >> eps;



auto start = chrono::high_resolution_clock::now();

	// Initialisations
	methode->Initialisation(b,A,sol0,r,results,methode);


	// On sauvegarde la solution
	methode->SaveSolution(0, r);



	while (r.norm()>eps && k<=k_max)
	{
		methode->calcul_sol(r);  			 //Appel de la fonction solution
		methode->SaveSolution(k,r);  	 //Sauvegarde de la norme de r

		cout << "k=" << k << "\n"<< endl;
		cout << "=======================" << endl;
		cout << "norme" <<r.norm() << endl;
		k+=1;
	}

	if (k>k_max)
	{
		cout << "Tolérance non atteinte :" << r.norm() << endl;
	}

	auto finish = chrono::high_resolution_clock::now();

	double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
	cout << "Cela a pris "<< t << " microsecondes" << endl;


	delete methode;
	return 0;
}


