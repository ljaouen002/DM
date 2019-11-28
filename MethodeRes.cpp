#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>
#include "Dense"
#include "Sparse"

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes()
{}

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

// Initialisation du nom du fichier
void MethodeRes::InitializeFileName(const string file_name)
{
  _file_out.open(file_name);
}

// Initialisation de vos différentes variables
void MethodeRes::Initialisation(VectorXd b, MatrixXd A, VectorXd sol0, VectorXd r,string results, MethodeRes* methode)
{
  _A=A;
  _b=b;
  _sol0=sol0;
  _sol=sol0;
  _r=r;
  _methode=methode;

  if (results.size() > 0)
  {
    _methode->InitializeFileName(results);
  }
}

void MethodeRes::SaveSolution(const int nb_iterations) //reste à déterminer ce qu'on trace
{
  _file_out << nb_iterations;
  _file_out << " " << _r.norm();
  _file_out << std::endl;
}



//Méthode de Jacobi

Jacobi::Jacobi()
{
  _D.setZero(_D.rows(),_D.cols());
  _E.setZero(_E.rows(),_E.cols());
  _F.setZero(_F.rows(),_F.cols());
  _M.setZero(_M.rows(),_M.cols());
  _N.setZero(_N.rows(),_N.cols());

  // Définition des matrices à utiliser dans le cas de Jacobi
  _D=_A.diagonal();   // Diagonale de A

  //Création de M=D^-1
  for (int i=0 ; i<_D.rows() ; ++i)
  {
    _M(i,i)=1/_D(i,i);
  }

  //Création de E et F
  for (int i=0 ; i<_A.rows() ; ++i)
  {
    for (int j=0 ; j<_A.cols() ; ++j)
    {
      if (i<j)
      {
        _E(i,j)=-_A(i,j);     // Partie triangulaire supérieure de A
      }
      else if (i>j)
      {
        _F(i,j)=-_A(i,j);     // Partie triangulaire inférieure de A
      }
    }
  }
  _N=_E+_F;

  _r=_b-_A*_sol0;
}

void Jacobi::calcul_sol()
{
  _sol=_M*_N*_sol+_M*_b;
  _r=_b-_A*_sol;
}


//Méthode du GPO

GPO::GPO()
{
  _r=_b-_A*_sol0;
}

void GPO :: calcul_sol()
{
  VectorXd _z;
  double _alpha;

  _z = _A *_r;
  _alpha = _r.dot(_r)/_z.dot(_r);
  _sol= _sol + _alpha*_r;
  _r=_b-_A*_z;
}


//Méthode du résidu minimum

Residu::Residu()
{
  _r=_b-_A*_sol0;
}


void Residu::calcul_sol()
{
  VectorXd _z;
  double _alpha;

  _z = _A * _r;
  _alpha = _r.dot(_z)/_z.dot(_z);
  _sol= _sol + _alpha*_r;
  _r=_b-_A*_sol;
}



/*
//Méthode de GMRes

GMRes::GMRes(double beta, SparseVector<double> r)
{
  _r=r;
  _beta=beta;
}

void GMRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A)
{
  _r=b-A*_sol0;
  //_beta = _r.norm(); pas besoin si beta existe déjà dans le main
}

void GMRes::Arnoldi(SparseVector<double> v, int N, SparseMatrix<double> A, SparseMatrix<double> v_arno, SparseMatrix<double>  H)
{
  //Définition des tailles des matrices H et v_arno, qui deviendront Hm et Vm
  v_arno.resize(N,N);
  H.resize(N,N);

  //definition des vecteurs et matrices locaux utilisés:
  SparseVector<double> wj(N); // Produit matrice vecteur A*vj, économiser des opérations
  SparseVector<double> vi(N); // Vecteur vi, économiser des opérations
  SparseVector<double> zj(N);
  SparseVector<double> Sk(N); // Vecteur résultant de la somme

  //Normalisation de v1
  v_arno.col(1)= (1/v.norm())*v;

  for (int j=1 ; j < N  ; j++)
  {
    wj= A*v_arno.col(j);
    for (int i=1 ; i < j  ; i++)
    {
      vi=v_ano.col(i);
      H(i,j)= wj.dot(vi);
    }

    for (int k=1 ; k < N  ; k++)
    {
      Sk = Sk + H(k,j)*v_arno.col(j);
    }

    zj = wj - Sk;
    H(j+1,j)= zj.norm();

    if (H(j+1,j)=0)
    {
      break;
    }
    v_arno.col(j+1) = zj / H(j+1,j) ;
  }
}

void GMRes::calcul_sol(SparseVector<double> b, SparseMatrix<double> A)
{}
*/

#define _METHODE_RES_CPP
#endif
