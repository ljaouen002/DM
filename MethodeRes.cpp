#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes(): _solution(0)
{}

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

// Initialisation de vos différentes variables
void MethodeRes::Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVextor<double> sol0, MethodeRes* methode)
{
  _A=A;
  _b=b;
  _sol0 = sol0;
  _sol = sol0;
  _methode = methode;

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






Jacobi::Jacobi(Eigen::SparseVector<double> r,Eigen::SparseVector<double> sol, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> F)
{
  _r=r;
  _sol=sol;
  _M=_D;
  _N=_E+_F;
}

void Jacobi::Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> F, SparseVector<double> r)
{
  _r=b-_A*_sol0;
}

void Jacobi::calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, SparseVector<double> r)
{
  _sol=_M*_N*_sol+_M*b; //Changer _M par _M^-1
  _r=b-_A*_sol;
}





Residu :: Residu(double alpha, Eigen::SparseVector<double> r)
{_r=r;
  _alpha =alpha
}

void Residu :: Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, SparseVector<double> r);
{
    _r=b-_A*_sol0;
}

void Residu :: calcul_sol(Eigen::SparseVector<double> q, Eigen::SparseMatrix<double> A, double alpha, SparseVector<double> r, SparseVector<double> z)
{
  z = _A * _r;
  _alpha = _r.dot(z)/z.dot(z);
  _sol= _sol + _alpha*_r;
  _r=b-_A*_sol;
}


#define _METHODE_RES_CPP
#endif
