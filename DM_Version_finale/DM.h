//#ifndef _DM_H

#include <vector>
#include <string>
#include <cmath>
#include <ostream>
#include "Dense"
#include "Sparse"
#include <fstream>
#include <iostream>


class Solveur

{
protected:
  // Vecteur solution X, résidu r, vecteur b et solution initiale X0
  Eigen::VectorXd _X, _b,_r,_X0;
  // A matrice de matrice market ou définie dans la question 1, Id la matrice identitée et M la matrice GSS
  Eigen::SparseMatrix<double> _A,_Id,_M;
  // Taile de la matrice carrée
  int _N,_m;
  std::ofstream _file_out;
  //


public:
  // Constructeur par défaut
  Solveur();
  // Destructeur par défaut
  virtual ~Solveur();
  // Fonction d'initialisation pour les matrices de MatrixMarket
  void Initialisation(std::string file_name, bool symetrique);
  // Fonction d'initialisation pour la question 1
  void Initialisation(int n, double alpha);
  // Méthode qui varie selon le solveur choisi
  virtual void method(int kmax, double eps)=0;
  // Initialise D, E et F
  void InitializeDEF();

};
// Classe fille pour le solveur GSS
class GSS : public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour le solveur Résidu minimum
class Res_min: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour le solveur gradient conjugué
class Grad_conj: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour le solveur gradient à pas optimal
class Grad_opt: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour le solveur GMRes
class GMRes: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour le solveur GMRes préconditionné
class GMResPreconditionne: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour la première partie de la partie théorique
class TH: public Solveur
{
public:
  void method(int kmax, double eps);
};

// Classe fille pour la partie théorique préconditionnée
class THCOND: public Solveur
{
public:
  void method(int kmax, double eps);
};
