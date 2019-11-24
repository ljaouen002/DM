#ifndef _METHODE_RES_H

#include "Dense"
#include "Sparse"
#include <fstream>

class MethodeRes
{
  private:
    // Vecteur initial et vecteur solution
    Eigen::VectorXd _sol0, _sol;
    // Pointeur vers le système d'EDO
    MethodeRes* _solution;
  protected:
    Eigen::SparseMatrix<double> _A;
    Eigen::SparseVector<double> _b;

  public:
    // Constructeur par défaut
    MethodeRes();
    // Destructeur par défaut
    virtual ~MethodeRes();
    // Initialisation
    virtual void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, MethodeRes* methode);
    // Calcul de x
    virtual void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A) =0;
    // Remplissage de la solution
    void SaveSolution(const int nb_iterations);
};

// Classe fille publique de MethodeRes
class Jacobi : public MethodeRes
{
  private:
    Eigen::SparseVector<double> _r;
    Eigen::SparseVector<double> _sol; //Lisa: Pour moi sol est déjà déclaré dans la classe mère
    Eigen::SparseMatrix<double> _M;
    Eigen::SparseMatrix<double> _N;
    //Eigen::SparseMatrix<double> _D, _E, _F;

  public:
    //constructeur
    Jacobi(Eigen::SparseVector<double> r, Eigen::SparseVector<double> sol);
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> F);
    void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);
};



// Classe fille publique de MethodeRes
class Residu : public MethodeRes
{
  private:
    double _alpha; //coefficient de descente
    Eigen::SparseVector<double> _r; // direction de descente, comme méthode de gradient, r=d
    //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3

  public:
    //constructeur
    Residu(double alpha, Eigen::SparseVector<double> r);
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> F);
    void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);
};


#define _METHODE_RES_H
#endif
