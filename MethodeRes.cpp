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
void MethodeRes::Initialisation(Eigen::SparseVector b, Eigen::SparseMatrix A, Eigen::SparseVextor sol0, MethodeRes* methode)
{
  _A=A.
  _b=b;
  _sol0 = sol0;
  _sol = sol0;
  _methode = methode;
  if (results.size() > 0)
  {
    _methode->InitializeFileName(results);
  }
}

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
  _sol=_M*_N*_sol+M*b; //Changer _M par _M^-1
  _r=b-A*_sol;
}

void MethodeRes::SaveSolution(const int nb_iterations) //reste à déterminer ce qu'on trace
{
  _file_out << nb_iterations;
  _file_out << " " << _r.norm();
  _file_out << std::endl;
}


#define _METHODE_RES_CPP
#endif
