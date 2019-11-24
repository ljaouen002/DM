#ifndef _TIME_SCHEME_H

#include "OdeSystem.h"

class TimeScheme
{
  private:
    // Pas de temps
    double _dt;
    // Temps en cours
    double _t;
    // Vecteur initial et vecteur solution
    Eigen::VectorXd _sol0, _sol;
    // Pointeur vers le système d'EDO
    OdeSystem* _sys;

  public:
    // Constructeur par défaut
    TimeScheme();
    // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
    // alors le compilateur en génère un implicitement.
    ~TimeScheme();
    // Initialisation de vos différentes variables
    void Initialize(double t0, double dt, Eigen::VectorXd & sol0, std::string name_file, OdeSystem* sys);
    // Enregistre la solution un fichier
    void SaveSolution();
    // Une étape du schéma en temps
    void Advance();
    // Permet de récupérer _sol
    const Eigen::VectorXd & GetIterateSolution() const;
};

#define _TIME_SCHEME_H
#endif
