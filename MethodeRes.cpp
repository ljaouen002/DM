#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes()
{
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  EIgen::MatrixXd I;
  double alpha(0.1);

  A.resize(10); B.resize(10); I.resize(10);
  
  for (int i=0 ; i<B.size() ; i++)
  {
    for (int j=0 ; j<B.size() ; j++)
    {
      I(i,i)=1.;
      B(i,j)=rand();
    }
  }
  A=alpha*I+B.transpose()*B;
}

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

MethodeRes::void calcul_residu(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::MatrixXd A)
{
  double residu;
  residu = b-A*x0 ;
}

Jacobi::Jacobi(double r, Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VextorXd x0)
{
  _r=r;
  _A=A;
  _b=b;
  _x0=x0;
}

void Jacobi::Initialisation()
{
  _r=calcul_residu(_x0,_b,_A)
}

#define _METHODE_RES_CPP
#endif