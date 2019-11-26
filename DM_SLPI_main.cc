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
SparseMatrix<double> Id(N,N), A(N,N), D(N,N), E(N,N), F(N,N), B(N,N);
SparseVector<double> sol(N), sol0(N);
double alpha(0.1);
SparseVector<double> r;
int userchoicemethode;

// Définition des matrices à utiliser globalement
A.setZero();

Id = MatrixXd::Identity(10,10);    //Matrice Identité
B=MatrixXd::Random(10,10);        //Matrice random B
//Lisa: Ce que je comprends du sujet, c'est que la matrice  ne doit contenir que des valeurs égale à 0 ou à 1, pas entre les deux.
A=alpha*Id+B.transpose()*B;        //Matrice A


cout << "------------------------------------" << endl;
cout << "Choississez la méthode de résolution : " << endl;
cout << "1) Jacobi"<< endl;
cout << "2) Gradient Pas Optimal" << endl;
cout << "3) Résidu Minimum" << endl;
cout << "4) GMRes" << endl;
//cin >> userchoicemethode;
userchoicemethode=1;
MethodeRes* methode(0);

switch(userchoicemethode)
{

	case 1: //Jacobi
	// Définition d'un pointeur de MethodeRes

	MethodeRes* methode = new Jacobi(r, sol);
	//Lisa: Pour moi il y a un problème de déclaration dans les variables de ta classe. Tu ne peux pas utiliser E, D et F dans ton constrcteur et dans ta classe sans le déclarer comme variable ptroctected. Donc pour moi, il faudrait mettre la ligne 52 à 67 dans la fonction initialized

	//Décalration des variables à l'intérieure du constrcteur, DONNER VALEUR

	// Nom du fichier solution
	string results = "solution_Jacobi.txt";

	E.setZero();
	F.setZero();

	// Définition des matrices à utiliser dans le cas de Jacobi
	D=diag(A);   // Diagonale de A, ne compire pas, retrouver fonction appropriée pour diagonale
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
	break;


	case 2: //GPO
			double alpha;
	MethodeRes* methode = new GPO(alpha, r);


			string results = "solution_GPO.txt";

	break;


	case 3: //Résidu

	double alpha;
	MethodeRes* methode = new Residu(alpha, r);

	string results = "solution_Residu.txt";
	break;


	default:
	cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
	exit(0);
}

// Initialisations
methode->Initialisation(b,A,sol0,methode);

// On sauvegarde la solution
methode->SaveSolution(nb_iterations);

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
