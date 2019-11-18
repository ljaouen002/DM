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
  double alpha(0.1);

  A.resize(10); B.resize(10);
  B=
}

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

MethodeRes::void calcul_residu(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::MatrixXd A)
{
  double residu;
  residu = b-A*x0 ;
}


#define _METHODE_RES_CPP
#endif
