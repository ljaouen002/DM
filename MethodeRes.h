#ifndef _METHODE_RES_H

#include "Dense"
#include <fstream>

class MethodeRes
{
  private:
    // Vecteur initial et vecteur solution
    Eigen::VectorXd _sol0, _sol;
    // Pointeur vers le système d'EDO
    MethodeRes* _solution;
  public:
    // Constructeur par défaut
    MethodeRes();
    // Destructeur par défaut
    virtual ~MethodeRes();
    // Initialisation
    virtual void Initialisation(Eigen::SparseVector b, Eigen::SparseMatrix A, Eigen::SparseVextor sol0, MethodeRes* methode);
    // Calcul de x
    virtual void calcul_sol(Eigen::SparseVector b, Eigen::SparseMatrix A) =0;
    // Remplissage de la solution
    void SaveSolution();
};

// Classe fille publique de MethodeRes
class Jacobi : public MethodeRes
{
  private:
    Eigen::SparseVector _r;
    Eigen::SparseVector _sol;
    Eigen::SparseMatrix _M;
    Eigen::SparseMatrix _N;
  public:
    Jacobi(Eigen::SparseVector<double> r, Eigen::SparseVector<double> sol);
    void Initialisation(Eigen::SparseVector b, Eigen::SparseMatrix A, Eigen::SparseVextor sol0, Eigen::SparseMatrix D, Eigen::SparseMatrix E, Eigen::SparseMatrix F);
    void calcul_sol(Eigen::SparseVector b, Eigen::SparseMatrix A);
};

#define _METHODE_RES_H
#endif
