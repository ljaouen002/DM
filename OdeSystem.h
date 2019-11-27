#ifndef _ODE_SYSTEM_H

#include "Dense"
#include <fstream>

class OdeSystem
{
  private:
    // Votre fonction f
    Eigen::VectorXd _f;
    // Écriture du fichier
    std::ofstream _file_out;
    
  public:
    // Constructeur par défaut
    OdeSystem();
    // Destructeur par défaut
    ~OdeSystem();
    // Initialiser le nom du fichier solution
    void InitializeFileName(const std::string file_name);
    // Sauvegarde la solution
    void SaveSolution(const double t, const Eigen::VectorXd & sol);
    // Pour récupérer _f
    const Eigen::VectorXd &  GetF()  const;
    // Pour construire _f en fonction de votre système
    void BuildF(const double t, const Eigen::VectorXd & sol);
};

#define _ODE_SYSTEM_H
#endif
