#ifndef _METHODE_RES_H

#include <fstream>
#include "Dense"
#include "Sparse"

class MethodeRes
{
protected:
  Eigen::SparseMatrix<double> _A;
  // Vecteur initial et vecteur solution
  Eigen::SparseVector<double> _sol0, _sol, _b, _r;
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
  virtual void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode)=0;
  // Calcul de x
  virtual void calcul_sol(Eigen::SparseVector<double>& _r) =0;
  // Remplissage de la solution
  void SaveSolution(const int nb_iterations, Eigen::SparseVector<double>& _r);


  //    void InitialisationMat(int& N,string name_Matrix, SparseMatrix<double>& A);
  // void InitialisationMat(std::string name_Matrix, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> b, Eigen::SparseVector<double> sol0, Eigen::SparseVector<double>& _r);

};

// Classe fille publique de MethodeRes
class Jacobi : public MethodeRes
{
private:
  Eigen::SparseMatrix<double> _M, _N;
public:
  Jacobi( Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> N);
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};

// Classe fille publique de MethodeRes
class GPO : public MethodeRes
{
private:
  //  double _alpha; //coefficient de descente
  //Eigen::VectorXd _z;
  //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3
public:
  GPO();
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};

// Classe fille publique de MethodeRes
class Residu : public MethodeRes
{
private:
  //Eigen::SparseMatrix<double> _M; //Matrice préconditionnement, dans question 3
public:
  Residu();
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};


// Classe fille publique de MethodeRes
class GMRes : public MethodeRes
{
private:
  Eigen::SparseMatrix<double> _Vm;
  Eigen::SparseMatrix<double> _Hm;
  //Eigen::SparseVector<double> _y;
  int _m;
public:
  //constructeur
  GMRes(Eigen::SparseMatrix<double> v_arno, Eigen::SparseMatrix<double> H, int m);
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void Arnoldi(Eigen::SparseVector<double> v, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> v_arno, Eigen::SparseMatrix<double> H);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};



// Classe fille publique de MethodeRes

class Residu_Precondi_gauche : public MethodeRes


{
private:
  Eigen::SparseMatrix<double> _M;
  Eigen::SparseMatrix<double> _F;
  Eigen::SparseMatrix<double> _E;
  Eigen::SparseMatrix<double> _D;
  Eigen::SparseMatrix<double> _D1;
  Eigen::SparseVector<double> _q;
  int _precondi;
public:

  Residu_Precondi_gauche(Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> F, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> D1, Eigen::SparseVector<double> q , int precondi);
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};


// Classe fille publique de MethodeRes
class Residu_Precondi_droite : public MethodeRes
{
private:
  Eigen::SparseMatrix<double> _M;
  Eigen::SparseMatrix<double> _F;
  Eigen::SparseMatrix<double> _E;
  Eigen::SparseMatrix<double> _D;
  Eigen::SparseMatrix<double> _D1;
  int _precondi;
public:
  Residu_Precondi_droite(Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> F, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> D1, int precondi);
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};


class Residu_Precondi_auto : public MethodeRes
{
private:

public:
  Residu_Precondi_auto();
  void Initialisation(Eigen::SparseVector<double> b, Eigen::SparseMatrix<double> A, Eigen::SparseVector<double> sol0 , Eigen::SparseVector<double>& _r, std::string results, MethodeRes* methode);
  void calcul_sol(Eigen::SparseVector<double>& _r);
};



#define _METHODE_RES_H
#endif

