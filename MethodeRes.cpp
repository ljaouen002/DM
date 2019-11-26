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
  _D=A.diagonal();   // Diagonale de A

  //Création de M=D^-1
  for (int i=0, i<D.rows() ; i++)
  {
    _M=1/-D(i,i);
  }

  for (int i=0 ; i<A.rows() ; i++)
  {
    for (int j=0 ; j<A.cols() ; j++)
    {
      if (i<j)
      {
        _E=-A(i,j);     // Partie triangulaire supérieure de A
      }
      else if (i>j)
      {
        _F=-A(i,j);     // Partie triangulaire inférieure de A
      }
    }
  }
  _N=_E+_F;

  _r=b-_A*_sol0;
}

void Jacobi::calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A)
{
  _sol=_M*_N*_sol+_M*b; //Changer _M par _M^-1
  _r=b-_A*_sol;
}


//Méthode du résidu minimum


Residu :: Residu(double alpha, Eigen::SparseVector<double> r)
{
  _r=r;
  _alpha =alpha;
}

void Residu :: Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, SparseVector<double> r)
{
  _r=b-_A*_sol0;
}


void Residu :: calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, /*double alpha, SparseVector<double> r,*/ SparseVector<double> z)
{
  z = _A * _r;
  _alpha = _r.dot(z)/z.dot(z);
  _sol= _sol + _alpha*_r;
  _r=b-_A*_sol;
}


GPO :: GPO(double alpha, Eigen::SparseVector<double> r)
{
  _r=r;
  _alpha =alpha;
}

void GPO :: Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, SparseVector<double> r)
{
  _r=b-_A*_sol0;
}


void GPO :: calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, /*double alpha, SparseVector<double> r,*/ SparseVector<double> z)
{
  z = _A * _r;
  _alpha = _r.dot(_r)/z.dot(_r);
  _sol= _sol + _alpha*_r;
  _r=b-_A*z;
}

#define _METHODE_RES_CPP
#endif
