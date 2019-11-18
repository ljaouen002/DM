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

// Calcul du résidu
MethodeRes::void calcul_residu(Eigen::VectorXd x, Eigen::VectorXd b, Eigen::MatrixXd A)
{
  double residu;
  residu = b-A*x ;
}

Jacobi::Jacobi(Eigen::VectorXd r, Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VextorXd x0, Eigen::VextorXd x)
{
  _r=r;
  _A=A;
  _b=b;
  _x0=x0;
  _x=x;
}

void Jacobi::Initialisation()
{
  _r(0)=calcul_residu(_x0,_b,_A)
}

void Jacobi::calcul_x()
{
  double eps(0.0001);
  int k(0); k_max(100);
  while (norme(_r)>eps && k<=k_max)
  {
    k+=1
    _x(k)=calcul_xk(Eigen::MatrixXd _A, Eigen::VectorXd _b)
    _r(k)=calcul_residu(Eigen::VectorXd _x, Eigen::VectorXd _b, Eigen::MatrixXd _A);
  }
  cout << _x << endl;
  if (k>k_max)
  {
    cout << "Tolérance non atteinte :" << norme(_r) << endl;
  }
}

#define _METHODE_RES_CPP
#endif
