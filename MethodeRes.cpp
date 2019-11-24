#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes()
{}

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

Jacobi::Jacobi(Eigen::SparseVector r,Eigen::SparseVextor sol, Eigen::SparseMatrix D, Eigen::SparseMatrix E, Eigen::SparseMatrix F)
{
  _r=r;
  _sol=sol;
  _M=D;
  _N=E+F;
}

void Jacobi::Initialisation(Eigen::SparseVector b, Eigen::SparseMatrix A, Eigen::SparseVector sol0, Eigen::SparseMatrix E, Eigen::SparseMatrix F)
{
  _r=b-A*sol0;
}

void Jacobi::calcul_sol(Eigen::SparseVector b, Eigen::SparseMatrix A)
{
  double eps(0.0001);
  int k(0); k_max(100);
  while (_r.norm()>eps || k<=k_max)
  {
    k+=1
    _sol=M*_N*_sol+M*b
    _r=b-A*_sol;
  }
  //cout << _sol << endl;
  if (k>k_max)
  {
    cout << "Tolérance non atteinte :" << _r.norm() << endl;
  }
}

#define _METHODE_RES_CPP
#endif
