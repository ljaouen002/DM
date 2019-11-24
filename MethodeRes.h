#ifndef _METHODE_RES_H

#include "Dense"
#include <fstream>

class MethodeRes
{
  public:
    // Constructeur par défaut
    MethodeRes();
    // Destructeur par défaut
    virtual ~MethodeRes();
    // Initialisation
    virtual void Initialisation(Eigen::SparseVector b, Eigen::SparseMatrix A, Eigen::SparseVextor sol0) =0;
    // Calcul de x
    virtual void calcul_sol() =0;
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
