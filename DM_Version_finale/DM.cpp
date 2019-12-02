//#ifndef _DM_CPP

#include <vector>
#include <string>
#include <cmath>
#include <ostream>
#include "Dense"
#include "Sparse"
#include "DM.h"
#include <fstream>
#include <iostream>
#include <iomanip>


using namespace Eigen;
using namespace std;


// Constructeur
Solveur::Solveur(){}
//Destructeur
Solveur::~Solveur(){}


void Solveur::Initialisation(string file_name,  bool symetrique)
{
  // Taille de la matrice et nombre d'éléments non nuls
  double taille,nb_Non_Nuls;
  // Entier pour le nombre d'éléments non nuls
  int non_Nuls;
  _m=5;
  // Chaine de caractère pour lire l'indice de la ligne, de la colonne et la valeur
  string l, c, val;
  vector<Triplet<double>> triplets;
  // Ouvre le fichier file_name.mtx
  ifstream mon_flux(file_name+".mtx");
  // Récupère la taille de la matrice et le nombre d'éléments non nuls
  mon_flux >> taille >> taille >> nb_Non_Nuls;
  cout << taille << endl;
  // converti la taille de la matrice en entier
  _N=int(taille);
  // redimensionne tous les vecteurs et les matrices utilisés
  _b.resize(_N);
  _X.resize(_N);
  _X0.resize(_N);
  _r.resize(_N);
  _A.resize(_N,_N);
  _Id.resize(_N,_N);
  // converti le nombre d'éléments non nuls en entier
  non_Nuls=int(nb_Non_Nuls);

  // Récupère les triplets pour remplir la matrice creuse A
  for (int k = 0; k < non_Nuls; k++)
  {
    mon_flux >> l;
    mon_flux >> c;
    mon_flux >> val;
    triplets.push_back({stoi(l)-1,stoi(c)-1,stod(val)});
    if (c != l && symetrique)
    {
      triplets.push_back({stoi(c)-1,stoi(l)-1,stod(val)});
    }
  }
  mon_flux.close();
  // Rempli la matrice A
  _A.setFromTriplets(triplets.begin(), triplets.end());
  // défini la matrice identité
  _Id.setIdentity();
  // défini le vecteur b
  for (int i = 0; i < _N; i++)
  {
    _b[i]=1.;
  }
  // défini la solution initiale
  for (int i = 0; i < _X0.size(); i++)
  {
    _X0[i]=1.;
  }
  _X=_X0;
  // défini le résidu initial
  _r=_b-_A*_X;
}



void Solveur::Initialisation(int n, double alpha)
{
  _N=n;
  _m=5;
  MatrixXd C(_N,_N), Unit(_N,_N);
  // dimensionne les vecteurs et matrices utilisés
  _b.resize(_N);
  _r.resize(_N);
  _X.resize(_N);
  _X0.resize(_N);
  _A.resize(_N,_N);
  _Id.resize(_N,_N);
  SparseMatrix<double> B;


  // C est une matrice de réels aléatoires entre -1 et 1
  C=MatrixXd::Random(_N,_N);
  // Unit est la matrice dont les coefficients sont tous égaux à 1
  Unit=MatrixXd::Constant(_N,_N,1.);
  // La matrice C est désormais une matrice de réels aléatoires entre 0 et 1
  C=1./2.*(Unit+C);
  // B est la représentation creuse de C
  B=C.sparseView();
  // définition de la matrice identité
  _Id.setIdentity();
  // Initialisation du vecteur b
  for (int i = 0; i < _N; i++)
  {
    _b[i]=1.;
  }
  // Initialisation de A
  _A= alpha*_Id+B;
  // Initialisation de la solution initiale
  for (int i = 0; i < _X0.size(); i++)
  {
    _X0[i]=1.;
  }
  _X=_X0;
  // Initialisation du résidu initiale
  _r=_b-_A*_X;
}



void Solveur::InitializeDEF()
{
  SparseMatrix<double> D1,D,E,F;
  D1.resize(_N,_N);
  D.resize(_N,_N);
  E.resize(_N,_N);
  F.resize(_N,_N);
  _M.resize(_N,_N);

  // D = Diaginale de A
  for (int i=0; i<_N; i++)
  {
    D.insert(i,i)=_A.coeffRef(i,i);
  }
  // F= - triangle inférieur de A
  F = - _A.triangularView<StrictlyUpper>();
  // F= - triangle supérieur de A
  E = - _A.triangularView<StrictlyLower>();

  // D1 est la matrice inverse de D
  for (int i = 0; i < _N; i++)
  {
      D1.insert(i,i)=1./D.coeff(i,i);
  }

  // définition de la matrice _M
  _M=(D-E)*D1*(D-F);
}



// Solveur GSS
void GSS::method (int kmax, double eps)
{

VectorXd y,v,z ;
double beta;
int k;

k=1;
beta = _r.norm();

// définition du solveur de la méthode LU
SparseLU<SparseMatrix<double>> solver;
solver.compute(_M);
//résout My=b
y=solver.solve(_b);SparseMatrix<double> D, SparseMatrix<double> F, SparseMatrix<double> E,

// Création et ouverture d'un fichier GSS.txt
_file_out.open("GSS"+to_string(_N)+".txt");
// écrit le résidu initial dans le fichier
_file_out<< k << "  " << beta << endl;
// Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
while ((beta>eps) && (k<kmax+1))
  {
    // z=AX
    z=_A*_X;
    // résout Mv=A*X
    v=solver.solve(z);
    // X(k+1) = X(k)-M^(-1)AX+M^(-1)b
    _X=_X-v+y;
    // Résidu à l'itération k
    _r=_b-_A*_X;
    // beta prend la valeur de la norme 2 du résidu
    beta = _r.norm();
    // écrit la valeur de k et la norme du résidu dans le fichier "GSS.txt"
    k+=1;
    _file_out<< k << "  " << _r.norm() << endl;
  }
  _file_out.close();
if (k>kmax)
  {
    cout << "Tolérance non atteinte" << beta << endl;
  }
}



// Solveur Residu minimum
void Res_min::method(int kmax, double eps)
{
  VectorXd z;
  double beta,alpha;
  int k;
  beta = _r.norm();
  k=1;
  // Création et ouverture d'un fichier solution Res_min.txt
  _file_out.open("Res_min"+to_string(_N)+".txt");
  // écrit le résidu initial dans le fichier
  _file_out<< k << "  " << beta << endl;
  // Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
  while ((beta>eps) && (k<kmax))
  {
    // z=Ar
    z=_A*_r;
    // Alpha vaut le rapport entre le ps de r et z et celui de z et z
    alpha= _r.dot(z)/(z.dot(z));
    // X(k+1)=X(k)+alpha*r
    _X=_X+alpha*_r;
    // r(k+1)=r(k)+alpha*z
    _r=_r-alpha*z;
    beta=_r.norm();
    // écrit la valeur de k et la norme du résidu dans le fichier "Res_min.txt"
    k=k+1;
    _file_out<< k << "  " << beta << endl;
  }
  _file_out.close();
  if (k>kmax)
  {
    cout << "Tolérance non atteinte" << _r.norm() << endl;
  }
}



// Solveur Gradient conjugué
void Grad_conj::method(int kmax, double eps)
{
  VectorXd p,z,rk;
  double beta, gamma,alpha;
  int k;
  p=_r;
  beta = _r.norm();
  k=1;

  // Création et ouverture du fichier
  _file_out.open("Grad_conj"+to_string(_N)+".txt");
  // écrit le résidu initial dans le fichier
  _file_out<< k << "  " << beta << endl;
  // Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
  while ((beta>eps) && (k<kmax))
  {
    // Conserve la valeur du résidu en entrée de boucle
    rk=_r;
    z=_A*p;
    // Redéfinition du pas de descente
    alpha= _r.dot(_r)/(z.dot(p));
    _X=_X+alpha*p;
    // calcul du nouveau résidu
    _r=_r-alpha*z;
    gamma=_r.dot(_r)/(rk.dot(rk));
    p=_r+gamma*p;
    beta=_r.norm();
    k +=1;
    // écrit la valeur de k et la norme du résidu dans le fichier "Grad_conj.txt"
    _file_out<< k << "  " << beta << endl;
  }
  _file_out.close();
  if (k>kmax)
  {
    cout << "Tolérance non atteinte" << _r.norm() << endl;
  }
}



void Grad_opt::method(int kmax, double eps)
{
  VectorXd z,rk;
  double beta, alpha;
  int k;
  beta = _r.norm();
  k=1;
  alpha=0.;
// Création et ouverture du fichier
  _file_out.open("Grad_opt"+to_string(_N)+".txt");
  _file_out<< k << "  " << beta << endl;
  // Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
  while ((beta>eps) && (k<kmax))
  {
    z=_A*_r;
    // calcul du nouveau pas de descente
    alpha= _r.dot(_r)/(z.dot(_r));
    _X=_X+alpha*_r;
    // calcul du nouveau résidu
    _r=_b-_A*_X;
    beta=_r.norm();
    k +=1;
    // écrit la valeur de k et la norme du résidu dans le fichier "Grad_opt.txt"
    _file_out<< k << "  " << beta << endl;
    cout << k << "  " << beta << endl;
  }
  _file_out.close();
  if (k>kmax)
  {
    cout << "Tolérance non atteinte" << _r.norm() << endl;
  }
}


// Solveur GMRes
void GMRes::method(int kmax, double eps)
{
  // Création des différentes matrices et vecteurs utiles dans cette méthode
  SparseMatrix<double> Hm,Qm,Rm;
  SparseMatrix<double> Vm;
  VectorXd y,e1,gm;
  SparseVector<double> zj,a,r;
  double beta1;
  int k;

  // Allocation des tailles des matrices et vecteurs
  r.resize(_N);
  zj.resize(_N);
  a.resize(_N);
  Vm.resize(_N,_m);
  Hm.resize(_m+1,_m);
  Qm.resize(_m+1,_m+1);
  gm.resize(_m+1);
  Rm.resize(_m+1,_m);
  y.resize(_m);
  e1.resize(_m+1);

  // Création du m-ième vecteur de la base canonique de taille m
  for (int i=1; i<e1.size(); i++)
  {
    e1.coeffRef(i)=0.;
  }
  e1.coeffRef(0)=1.;

  beta1=_r.norm();
  k=1;
  // Creation et ouverture du fichier GMRES.txt
  _file_out.open("GMRes"+to_string(_N)+".txt");
    _file_out<< k << "  " << beta1 << endl;

  // Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
  while ((beta1>eps) && (k<kmax+1))
  {
    // Réalise une itération d'Arnoldi
    r=_r.sparseView();
    Vm.col(0)=r/r.norm();
    for (int j = 0; j < _m-1; j++)
    {
      a.setZero();
      for (int i = 0; i < j+1; i++)
      {
        Hm.coeffRef(i,j)=(_A*Vm.col(j)).dot(Vm.col(i));
        a+=Hm.coeff(i,j)*Vm.col(i);
      }
      zj=(_A*Vm.col(j))-a;
      Hm.coeffRef(j+1,j)=zj.norm();
      if (Hm.coeff(j+1,j)<pow(10,-14))
      {
        Hm.coeffRef(j+1,j)=0;
        break;}
      Vm.col(j+1)=zj/Hm.coeff(j+1,j);
    }
    // Fin de Arnoldi

    // Décompostion QR de Hm
    SparseQR<SparseMatrix<double>,COLAMDOrdering<int>> solver;
    Hm.makeCompressed();
    solver.compute(Hm);
    Qm=solver.matrixQ();
    Rm=solver.matrixR();
    gm=beta1*Qm.transpose()*e1;
    solver.compute(Rm);
    // Résolution de Rm*y=gm
    y=solver.solve(gm);

    _X=_X+Vm*y;
    // calcul du nouveau résidu
    _r=_r-_A*Vm*y;

    beta1=_r.norm();
    k+=1;

    // écriture dans le fichier
    _file_out<< k << "  " << beta1 << endl;
  }
  _file_out.close();
  if (k>kmax)
    {
      cout << "Tolérance non atteinte" << _r.norm() << endl;
    }
}



void GMResPreconditionne::method(int kmax, double eps)
{
  // Création des différentes matrices et vecteurs utiles dans cette méthode
  SparseMatrix<double> Hm,Qm,Rm;
  SparseMatrix<double> Vm;
  VectorXd y,e1,gm,q,w,w1,l;
  SparseVector<double> zj,a,r,q1,w2,w3;
  double beta1,b2;
  int k;

  // Allocation des matrices et des vecteurs
  r.resize(_N);
  q.resize(_N);
  q1.resize(_N);
  zj.resize(_N);
  a.resize(_N);
  Vm.resize(_N,_m);
  Hm.resize(_m+1,_m);
  Qm.resize(_m+1,_m+1);
  gm.resize(_m+1);
  Rm.resize(_m+1,_m);
  y.resize(_m);
  e1.resize(_m+1);
  w.resize(_N);
  w1.resize(_N);
  w2.resize(_N);
  w3.resize(_N);
  l.resize(_N);

  // Création du m-ième vecteur de la base canonique de taille m
  for (int i=1; i<e1.size(); i++)
  {
    e1.coeffRef(i)=0.;
  }
  e1.coeffRef(0)=1.;

  beta1=_r.norm();
  // décomposition LU pour résoudre Mq=r
  SparseLU<SparseMatrix<double>> solver;
  solver.compute(_M);
  q=solver.solve(_r);
  // norme du "faux résidu"
  b2=q.norm();
  k=1;
  // Ouverture et création du fichier GMRes_PrecC.txt
  _file_out.open("GMRes_PreC"+to_string(_N)+".txt");
  _file_out<< k << "  " << beta1 << endl;
  while ((beta1>eps) && (k<kmax+1))
  {
    // Arnoldi pour la base (q,(M^-1A)q, ... ,(M^-1A)^(m-1)q)
    q1=q.sparseView();
    Vm.col(0)=q1/q1.norm();
    for (int j = 0; j < _m-1; j++)
    {
      a.setZero();
      w3=_A*Vm.col(j);
      w1=VectorXd(w3);
      w=solver.solve(w1);
      w2=w.sparseView();
      for (int i = 0; i < j+1; i++)
      {
        Hm.coeffRef(i,j)=(w2).dot(Vm.col(i));
        a+=Hm.coeff(i,j)*Vm.col(i);
      }
      zj=w2-a;
      Hm.coeffRef(j+1,j)=zj.norm();
      if (Hm.coeff(j+1,j)<pow(10,-14))
      {
        Hm.coeffRef(j+1,j)=0;
        break;}
      Vm.col(j+1)=zj/Hm.coeff(j+1,j);
    }
    // Fin de Arnoldi

    // Décomposition QR
    SparseQR<SparseMatrix<double>,COLAMDOrdering<int>> solver2;
    Hm.makeCompressed();
    solver2.compute(Hm);
    Qm=solver2.matrixQ();
    Rm=solver2.matrixR();
    gm=b2*Qm.transpose()*e1;
    SparseQR<SparseMatrix<double>,COLAMDOrdering<int>> solver3;
    solver3.compute(Rm);
    // Résout Rm*y=gm
    y=solver3.solve(gm);

    _X=_X+Vm*y;
    // Calcul du nouveau "vrai" résidu
    _r=_b-_A*_X;
    // Calcul du nouveau "faux" résidu
    q=solver.solve(_r);

    beta1=_r.norm();
    b2=q.norm();
    k+=1;
    // Ecrit dans le fichier
    _file_out<< k << "  " << beta1 << endl;
  }
  _file_out.close();
  if (k>kmax)
    {
      cout << "Tolérance non atteinte" << _r.norm() << endl;
    }
}

void TH::method(int kmax, double eps)
{
  // Initialisation des variables nécessaires
  SparseMatrix<double> TA;
  VectorXd rk,rk1,p,z,TAb,q,w,x;
  double Beta,Alpha,Gamma,Norm;
  int k;
  // Création d'un fichier Th.txt pour accueillir les résidus à chaque itération
  _file_out.open("Th"+to_string(_N)+".txt");
  // Vrai résidu
  q = _b- _A*_X;
  // Calcul de la transposée de A
  TA = _A.transpose();
  // Calcul second membre
  TAb = TA * _b;
  x = _A * _X;
  x = TA * x;
  // Résidu du système préconditionné
  rk = TAb - x;
  p = rk;
  // Beta: norme 2 du vrai résidu
  Beta = q.norm();
  k = 0;

  while ((Beta>eps) && (k<kmax+1))
  {
    // w = Ap
    w = _A * p;
    // z1 = TA * A *p
    z = _A * p;
    z = TA *z;
    //Alpha = <rk,rk> / <z,p>
    Alpha = (rk.dot(rk))/(z.dot(p));
    // r(k+1) = r(k) - Alpha*z1
    rk1 = rk - Alpha * z; // Calcul de r(k+1)
    // Gamma = <r(k+1),r(k+1)> / <r(k),r(k)>
    Gamma = (rk1.dot(rk1))/(rk.dot(rk));
    // p = r(k+1) + Gamma * p
    p = rk1 + Gamma*p;
    // q =  q - Alpha * w
    q = q - Alpha * w;
    //Beta = norme du vrai résidu
    Beta = q.norm();
    // x = x - Alpha*p
    _X = _X - Alpha * p;
    k +=1;
    _file_out<< k << "  " << Beta << endl;
    //mise à jour de r(k)
    rk = rk1;
  }
  _file_out.close();
}

void THCOND::method(int kmax, double eps)
{
  // définition des variables nécessaires au calcul
  SparseMatrix<double> TA,TM;
  VectorXd rk,rk1,p,z1,z2,TMAb,y,z,TT;
  double Beta,Beta2, Alpha, Gamma, ratio;
  int k;
  _file_out.open("ThCond"+to_string(_N)+".txt");
  // Dimensionnent des matrices utilisées
  TM.resize(_N,_N);
  // Transposée de A
  TA = _A.transpose();
  // Initialisation k
  k=0;


  // Transposée de M
  TM = _M.transpose();

  //Définition des solveurs
  SparseLU<SparseMatrix<double>> solver,solver1;
  solver.compute(_M);
  solver1.compute(TM);

  // y = M*X
  y = _M* _X;

  // z: M(-T)A(T)AM(-1)Y
  z = _M * y;
  z = _A * z;
  z = TA * z;

  //Définition source terme
  TMAb = TA * _b;
  // TMAb = M^(-T) A^(T) b
  TT = solver1.solve(TMAb);
  TMAb = TT;

  //résout z = M(-1) * z
  z1=solver1.solve(y);
  z = z1;
  // résidu préconditionné
  rk = z - TMAb;
  // p = r(0)
  p = rk;
  // Beta: Norme 2 du résidu
  Beta = rk.norm();
  Beta2 = (_b - _A*_X).norm();
  // ratio: on divise les normes des résidus
  // du système de base et du système préconditionné
  ratio = Beta/Beta2;
  _file_out<< 0 << "  " << Beta << endl;

  while ((Beta2>eps*ratio) && (k<kmax+1))
  {

    //z1 = M(-T)A(T)AM(-1)p
    z1 = solver.solve(p);
    z1 = _A * z1;
    z1 = TA * z1;
    z2 = solver1.solve(z1);
    z1 = z2;

    // Alpha = <rk,rk>/<z1,p>
    Alpha = (rk.dot(rk))/(z1.dot(p));
    //r(k+1) = r(k) - Alpha * z1;
    rk1 = rk - Alpha * z1;
    //Gamma = <r(k+1),r(k+1)>/<rk,rk>
    Gamma = (rk1.dot(rk1))/(rk.dot(rk));

    // p = r(k+1) + Gamma*p
    p = rk1 + Gamma*p;

    //Beta : norme du résidu
    Beta = rk.norm();
    Beta2 = Beta / ratio;
    // y(k+1) = y(k) - Alpha*p
    y = y - Alpha * p;
    k +=1;
    _file_out<< k << "  " << Beta2 << endl;
    // mise à jour de r(k)
    rk = rk1;
  }
  _file_out.close();
}



//  #define _DM_CPP
//  #endif
