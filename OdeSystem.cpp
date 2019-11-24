#ifndef _ODE_SYSTEM_CPP

#include "OdeSystem.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
OdeSystem::OdeSystem()
{}

// Destructeur par défaut
OdeSystem::~OdeSystem()
{}

// Initialisation du nom du fichier
void OdeSystem::InitializeFileName(const std::string file_name)
{
  _file_out.open(file_name);
}

// Renvoie le vecteur _f
const VectorXd & OdeSystem::GetF() const
{
  return _f;
}


// Construit le vecteur f(t, sol)
void OdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = sol;
}

// Enregistre la solution
// Pour le moment : sol_1, sol_2 ...
void OdeSystem::SaveSolution(const double t, const VectorXd & sol)
{
  _file_out << t;
  for (int i = 0 ; i < sol.rows() ; ++i)
  {
    _file_out << " " << sol(i);
  }
  _file_out << std::endl;
}


#define _ODE_SYSTEM_CPP
#endif
