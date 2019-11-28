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
void MethodeRes::InitializeFileName(const std::string file_name)
{
  _file_out.open(file_name);
}

// Initialisation de vos différentes variables
void MethodeRes::Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, std::string results, MethodeRes* methode)
{
  _A=A;
  _b=b;
  _sol0=sol0;
  _sol=sol0;
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

Jacobi::Jacobi(Eigen::SparseVector<double> r,Eigen::SparseVector<double> sol)
{
  _r=r;
  _sol=sol;
}

void Jacobi::Initialisation(int N, Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0)
{
  _D.setZero(_D.rows(),_D.cols());
  _E.setZero(_E.rows(),_E.cols());
  _F.setZero(_F.rows(),_F.cols());
  _M.setZero(_M.rows(),_M.cols());
  _N.setZero(_N.rows(),_N.cols());
  
  // Définition des matrices à utiliser dans le cas de Jacobi
  _D=_A.diagonal();   // Diagonale de A

  //Création de M=D^-1
  for (int i=0, i<_D.rows() ; i++)
  {
    _M=1/-D(i,i);
  }

  for (int i=0 ; i<_A.rows() ; i++)
  {
    for (int j=0 ; j<_A.cols() ; j++)
    {
      if (i<j)
      {
        _E=-_A(i,j);     // Partie triangulaire supérieure de A
      }
      else if (i>j)
      {
        _F=-_A(i,j);     // Partie triangulaire inférieure de A
      }
    }
  }
  _N=_E+_F;

  _r=_b-_A*_sol0
}

void Jacobi::calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A)
{
  _sol=_M*_N*_sol+_M*_b; //Changer _M par _M^-1
  _r=_b-_A*_sol;
}


//Méthode du résidu minimum

Residu :: Residu(double alpha, Eigen::SparseVector<double> r)
{
  _r=r;
  _alpha =alpha;
}

void Residu :: Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, SparseVector<double> r)
{
  _r=_b-_A*_sol0;
}


void Residu :: calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A /*double alpha, SparseVector<double> r,*/ )
{
  SparseVector<double> z;
  z = _A * _r;
  _alpha = _r.dot(z)/z.dot(z);
  _sol= _sol + _alpha*_r;
  _r=_b-_A*_sol;
}


//Méthode du GPO

GPO :: GPO(double alpha, Eigen::SparseVector<double> r)
{
  _r=r;
  _alpha =alpha;
}

void GPO :: Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, SparseVector<double> r)
{
  _r=b-_A*_sol0;
}


void GPO :: calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, /*double alpha, SparseVector<double> r,*/)
{
  SparseVector<double> z;
  z = _A * _r;
  _alpha = _r.dot(_r)/z.dot(_r);
  _sol= _sol + _alpha*_r;
  _r=_b-_A*z;
}


//Méthode de GMRes

GMRes :: GMRes(double beta, Eigen::SparseVector<double> r)
{
  _r=r;
  _beta=beta;
}

void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, double beta)
{
  _r=_b-_A*_sol0;
  _beta = _r.norm();
}

void Arnoldi(Eigen::SparseVector<double> v, int N, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> v_arno, SparseMatrix<double>  H)
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


void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A)
{}







#define _METHODE_RES_CPP
#endif
