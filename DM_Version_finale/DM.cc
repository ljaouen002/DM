
#include <vector>
#include <string>
#include <cmath>
#include <ostream>
#include "Dense"
#include "Sparse"
#include "DM.h"
#include <fstream>
#include <iostream>
using namespace std;
using namespace Eigen;


int main ()
{
  int n;
  Solveur* Solv;
  string name_file;
  bool booleen;
  int userChoicee;
  int userChoice;
  int userChoiceee;
  int userChoiceeee;

    cout << "------------------------------------" << endl;
    cout << "Choississez le solveur : " << endl;
    cout << "1) GSS" << endl;
    cout << "2) Résidu minimum" << endl;
    cout << "3) Gradient conjugué" << endl;
    cout << "4) GMRes" << endl;
    cout << "5) GMRes Préconditionné" << endl;
    cout << "6) Gradient à pas optimal" << endl;
    cout << "7) Partie théorique 1" << endl;
    cout << "8) Partie théorique préconditionnée" << endl;
    cin >> userChoice;



    switch(userChoice)
  {
  case 1:
  Solv = new GSS();

  break;

  case 2:
  Solv = new Res_min();

  break;

  case 3:
  Solv = new Grad_conj();

  break;

  case 4:
  Solv = new GMRes();
  break;

  case 5:
  Solv = new GMResPreconditionne();

  break;
  case 6:
  Solv = new Grad_opt();

  case 7:
  Solv = new TH();
  break;

  case 8:
  Solv = new THCOND();


  break;

  default:
  cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
  exit(0);
  }


  cout << "----------------------------"  << endl;
  cout << "Quelle matrice voulez vous" << endl;
  cout << "1) Matrice An " << endl;
  cout << "2) Matrice du Matrix Market" << endl;
  cin >> userChoicee;

  switch(userChoicee)
  {


    case 1:
cout << "----------------------------------------" << endl;
    cout << "Donner la valeur de n:" << endl;
    cin >> n ;
    Solv->Initialisation(n,3*n);
    break;

    case 2:

cout << "----------------------------------------" << endl;

    cout << "Quelle matrice du matrixMarket" << endl;
    cout << "1) s3rmt3m3 qui est SDP, K=2.6*10^(10)" << endl;
    cout << "2) bcsstm24 qui est SDP, K=1.8*10^(13)" << endl;
    cout << "3) af23560 qui est non sym, K non donné" << endl;
    cout << "4) mcca qui est non sym, K=3.6*10^(17)" << endl;
cin >> userChoiceee;

  switch(userChoiceee)

  { case 1:
    name_file="s3rmt3m3";
  break;

  case 2:
  name_file="bcsstm24";
  break;

  case 3:
  name_file="af23560";
  break;

  case 4:
  name_file="mcca";
  break;

  default:
  cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
  exit(0);
}
cout << "----------------------------------------" << endl;
 cout << "Cette matrice est elle symétrique?" << endl;
 cout << "1)Oui" << endl;
 cout << "2)Non " << endl;
 cin >> userChoiceeee;

 switch(userChoiceeee)

 {
   case 1:
   booleen=true;
   Solv->Initialisation(name_file,booleen);
   break;

   case 2:
   booleen=false;
   Solv->Initialisation(name_file,booleen);
   break;

  default:
  cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
  exit(0);

        }

}

int nbite;
cout << "----------------------------------------" << endl;
 cout << "Combien d'itérations voulez vous?" << endl;
 cin >> nbite;

int tol;
cout <<"------------------------------------------"<< endl;
cout << "quelle puissance pour la tolérance (valeur à rentrer avec le -)" << endl;
cin >> tol;
Solv->InitializeDEF();
Solv->method(nbite, pow(10,tol));
delete Solv;
return 0;
}
