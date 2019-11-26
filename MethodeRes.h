#ifndef _METHODE_RES_H

#include <fstream>
#include "Dense"
#include "Sparse"

class MethodeRes
{
  protected:
    Eigen::SparseMatrix<double> _A;
    Eigen::SparseVector<double> _b;
    // Écriture du fichier
    std::ofstream _file_out;
    // Vecteur initial et vecteur solution
    Eigen::VectorXd _sol0, _sol;
    MethodeRes* _methode;
  public:
    // Constructeur par défaut
    MethodeRes();
    // Destructeur par défaut
    virtual ~MethodeRes();
    // Initialiser le nom du fichier solution
    void InitializeFileName(const std::string file_name);
    // Initialisation
    virtual void Initialisation(int N, Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0, std::string results, MethodeRes* methode);
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
    //Eigen::SparseVector<double> _sol; //Lisa: Pour moi sol est déjà déclaré dans la classe mère
    Eigen::SparseMatrix<double> _M, _N, _D, _E, _F;

  public:
    //constructeur
    Jacobi(Eigen::SparseVector<double> r, Eigen::SparseVector<double> sol);
    void Initialisation(int N, Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0);
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
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0);
    void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);
};

// Classe fille publique de MethodeRes
class GPO : public MethodeRes
{
  private:
    double _alpha; //coefficient de descente
    Eigen::SparseVector<double> _r; // direction de descente, comme méthode de gradient, r=d
    //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3

  public:
    //constructeur
    GPO(double alpha, Eigen::SparseVector<double> r);
    void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0);
    void calcul_sol(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A);
};

/*
// Classe fille publique d’OdeSystem
class SecondExampleOdeSystem : public OdeSystem
{
  public:
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = (-y,x)
};

// Classe fille publique d’OdeSystem
class ThirdExampleOdeSystem : public OdeSystem
{
  public:
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = tX²
};

// Classe fille publique d’OdeSystem
class LotkaVolterraOdeSystem : public OdeSystem
{
  private:
    //Paramètre Lotka Volterra
    double _a;
    double _b;
    double _c;
    double _d;
  public:
    LotkaVolterraOdeSystem(double a, double b, double c, double d);
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = (x(a-by),y(cx-d))
};

// Classe fille publique d’OdeSystem
class PendulumOdeSystem : public OdeSystem
{
  private:
    //Paramètre Lotka Volterra
    double _l;
    double _m;
    double _k;
  public:
    PendulumOdeSystem(double l, double m); //_k=0
    PendulumOdeSystem(double l, double m, double k);
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = (y,-g/l*sin(x)-k/ml²*y)
    void SaveSolution(const double t, const Eigen::VectorXd & sol);
};

// Classe fille publique d’OdeSystem
class CircuitRLCOdeSystem : public OdeSystem
{
  private:
    double _R;
    double _C;
    double _q0;
  public:
    CircuitRLCOdeSystem(double R, double C, double q0);
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = -(1/RC)*X
    void SaveSolution(const double t, const Eigen::VectorXd & sol);
};
*/
#define _METHODE_RES_H
#endif
