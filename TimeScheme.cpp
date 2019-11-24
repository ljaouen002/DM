#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _sys(0)
{}

// Destructeur
TimeScheme::~TimeScheme()
{}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(double t0, double dt, VectorXd & sol0, string results, OdeSystem* sys)
{
  _dt = dt;
  _t = t0 ;
  _sol0 = sol0;
  _sol = sol0;
  _sys = sys;
  if (results.size() > 0)
  {
    _sys->InitializeFileName(results);
  }
}

// Schéma en temps par défaut : ici Euler Explicite
// Avancer d'un pas de temps
void TimeScheme::Advance()
{
  _sys->BuildF(_t, _sol);
  _sol += _dt*_sys->GetF();
  _t += _dt;
}

// Enregistre la solution : fait appel à OdeSystem car la solution
// que l'on souhaite sauvegarder peut être différente de _sol SaveSolution
// le système
void TimeScheme::SaveSolution()
{
  _sys->SaveSolution(_t, _sol);
}

// Renvoie _sol (pratique pour calculer l'ordre de votre méthode)
const VectorXd & TimeScheme::GetIterateSolution() const
{
  return _sol;
}

#define _TIME_SCHEME_CPP
#endif
