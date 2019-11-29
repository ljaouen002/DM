#ifndef _METHODE_RES_H

#include <fstream>
#include "Dense"
#include "Sparse"

class MethodeRes
{
  protected:
    Eigen::SparseMatrix<double> _A, _b;
    // Vecteur initial et vecteur solution
    Eigen::SparseVector<double> _sol0, _sol, _r;
    // Écriture du fichier
    std::ofstream _file_out;
    MethodeRes* _methode;
  public:
    // Constructeur par défaut
    MethodeRes();
    // Destructeur par défaut
    virtual ~MethodeRes();
    // Initialiser le nom du fichier solution
    void InitializeFileName(const std::string file_name);
    // Initialisation
    virtual void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double> r, std::string results, MethodeRes* methode)=0;
    // Calcul de x
    virtual void calcul_sol() =0;
    // Remplissage de la solution
    void SaveSolution(const int nb_iterations);
};

// Classe fille publique de MethodeRes
class Jacobi : public MethodeRes
{
  private:
    Eigen::SparseMatrix<double> _D, _E, _F, _M, _N;
  public:
    Jacobi(Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> F, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> N);
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double> r, std::string results, MethodeRes* methode);
    void calcul_sol();
};

// Classe fille publique de MethodeRes
class GPO : public MethodeRes
{
  private:
    double _alpha; //coefficient de descente
    //Eigen::VectorXd _z;
    //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3
  public:
    GPO();
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double> r, std::string results, MethodeRes* methode);
    void calcul_sol();
};

// Classe fille publique de MethodeRes
class Residu : public MethodeRes
{
  private:
    //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3
  public:
<<<<<<< HEAD
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double> r, std::string results, MethodeRes* methode);
    void calcul_sol();
=======
    void Initialisation(Eigen::VectorXd b, Eigen::MatrixXd A, Eigen::VectorXd sol0 , Eigen::VectorXd r, std::string results, MethodeRes* methode);
void calcul_sol();
>>>>>>> 53845d8384c4d3c5a7a406de1396d39fa4844188
};


/*
// Classe fille publique de MethodeRes
class GMRes : public MethodeRes
{
  private:
    double _beta;
  public:
    //constructeur
    GMRes(double beta, Eigen::SparseVector<double> r);
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0);
    void Arnoldi(Eigen::SparseVector<double> v, int N, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> v_arno, SparseMatrix<double> H);
    void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);
};
*/
#define _METHODE_RES_H
#endif
