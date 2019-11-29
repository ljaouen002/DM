#ifndef _METHODE_RES_H

#include <fstream>
#include "Dense"
#include "Sparse"

class MethodeRes
{
  protected:
    Eigen::MatrixXd _A, _b;
    // Vecteur initial et vecteur solution
    Eigen::VectorXd _sol0, _sol, _r;
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
    virtual void Initialisation(Eigen::VectorXd b, Eigen::MatrixXd A, Eigen::VectorXd sol0 , Eigen::VectorXd r, std::string results, MethodeRes* methode)=0;
    // Calcul de x
    virtual void calcul_sol() =0;
    // Remplissage de la solution
    void SaveSolution(const int nb_iterations);

    Eigen::VectorXd Get_r(){return _r;};


};

// Classe fille publique de MethodeRes
class Jacobi : public MethodeRes
{
  private:
    Eigen::MatrixXd _D, _E, _F, _M, _N;
  public:
    Jacobi(Eigen::MatrixXd D, Eigen::MatrixXd F, Eigen::MatrixXd E, Eigen::MatrixXd M, Eigen::MatrixXd N);
    void Initialisation(Eigen::VectorXd b, Eigen::MatrixXd A, Eigen::VectorXd sol0 , Eigen::VectorXd r, std::string results, MethodeRes* methode);
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
    //constructeur
    GPO();
    void Initialisation(Eigen::VectorXd b, Eigen::MatrixXd A, Eigen::VectorXd sol0 , Eigen::VectorXd r, std::string results, MethodeRes* methode);
    void calcul_sol();
};

// Classe fille publique de MethodeRes
class Residu : public MethodeRes
{
  private:
    double _alpha; //coefficient de descente
    Eigen::VectorXd _z;
    //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3
  public:
    Residu();
    void Initialisation(Eigen::VectorXd b, Eigen::MatrixXd A, Eigen::VectorXd sol0 , Eigen::VectorXd r, std::string results, MethodeRes* methode);
    void calcul_sol();
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
