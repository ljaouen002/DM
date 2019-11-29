#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>
#include "Dense"
#include "Sparse"

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes()
{}
//
// MethodeRes::MethodeRes(VectorXd b, MatrixXd A, VectorXd sol0, VectorXd r,string results, MethodeRes* methode)
//   {
//     _A=A;
//     _b=b;
//     _sol0=sol0;
//     _sol=sol0;
//     _r=r;
//     _methode=methode;
//
//     if (results.size() > 0)
//     {
//       _methode->InitializeFileName(results);
//     }
//   }

// Destructeur par défaut
MethodeRes::~MethodeRes()
{}

// Initialisation du nom du fichier
void MethodeRes::InitializeFileName(const string file_name)
{
  _file_out.open(file_name);
}

// Initialisation de vos différentes variables
void MethodeRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
  _A=A;
  _b=b;
  _sol0=sol0;
  _sol=sol0;
  //_r=r;
  _methode=methode;

  if (results.size() > 0)
  {
    _methode->InitializeFileName(results);
  }
}

void MethodeRes::SaveSolution(const int nb_iterations, Eigen::SparseVector<double>& _r) //reste à déterminer ce qu'on trace
{
  _file_out << nb_iterations;
  _file_out << " " << _r.norm();
  _file_out << std::endl;
}



//Méthode de Jacobi

Jacobi::Jacobi(SparseMatrix<double> D, SparseMatrix<double> F, SparseMatrix<double> E, SparseMatrix<double> M, SparseMatrix<double> N)
{
  _D=D;
  _F=F;
  _E=E;
  _M=M;
  _N=N;
}

void Jacobi :: Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)

{
  MethodeRes::Initialisation(b,A,sol0,_r,results,methode);
  // _A=A;
  // _b=b;
  // _sol0=sol0;
  // _sol=sol0;
  // _r=r;
  // _methode=methode;
  //
  // if (results.size() > 0)
  // {
  //   _methode->InitializeFileName(results);
  // }


  _D=0*_A;
  _E=0*_A;
  _F=0*_A;
  _M=0*_A;
  _N=0*_A;

  // Définition des matrices à utiliser dans le cas de Jacobi
//  _D=_A.diagonal();   // Diagonale de A


  //Création de M=D^-1
  for (int i=0 ; i<_D.rows() ; ++i)
  {
    //  cout << _D.row() << endl;
    _D.coeffRef(i,i) = _A.coeffRef(i,i);
    _M.coeffRef(i,i)=1/_D.coeffRef(i,i);
  }

  //Création de E et F
  for (int i=0 ; i<_A.rows() ; ++i)
  {
    for (int j=0 ; j<_A.cols() ; ++j)
    {
      if (i<j)
      {

        _E.coeffRef(i,j)=-_A.coeffRef(i,j);     // Partie triangulaire supérieure de A
      //  cout << "E"<<_E.coeffRef(i,j) << endl;
      }
      else if (i>j)
      {
        _F.coeffRef(i,j)=-_A.coeffRef(i,j);     // Partie triangulaire inférieure de A
      //  cout << "F" << _F.coeffRef(i,j) << endl;
      }
    }
  }
  _N=_E+_F;

  _r=_b-_A*_sol0;

  cout << "F" << _F << endl;
  cout << "A" << _A << endl;
  cout << "_r=" << _r << endl;


	//	cout <<_r.norm() << " " <<  "1" << endl;>>>>>>> 409d716ac4dbc893d33b79cd73a3ca18ae722669>>>>>>> 53845d8384c4d3c5a7a406de1396d39fa4844188
}


void Jacobi::calcul_sol(Eigen::SparseVector<double>& _r)
{
  _sol=_M*_N*_sol+_M*_b;

  _r=_b-_A*_sol;
//  cout <<_N.norm() << " " <<  "1" << endl;
//  cout <<_r.norm() << " " <<  "2" << endl;

}



//Méthode du GPO

GPO::GPO()
{

}


void GPO :: Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
    MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

    _r=_b-_A*_sol0;
}

void GPO :: calcul_sol(Eigen::SparseVector<double>& _r)
{
  SparseVector<double> z;
  //double _alpha;

  z = _A *_r;
  _alpha = _r.dot(_r)/z.dot(_r);
  _sol= _sol + _alpha*_r;
 _r=_r-_alpha*z;


}


//Méthode du résidu minimum


void Residu::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
  MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

  _r=_b-_A*_sol;
}

void Residu::calcul_sol(Eigen::SparseVector<double>& _r)
{
  double _alpha;

  SparseVector<double> _z;


  _z=_A*_r;
  _alpha=_r.dot(_z)/_z.dot(_z);
  _sol=_sol+_alpha*_r;
  _r=_r-_alpha*_z;
}



/*
//Méthode de GMRes
GMRes::GMRes(double beta, Eigen::SparseVector<double>& _r)
{
  _r=r;
  _beta=beta;
}
void GMRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A)
{
  _r=b-A*_sol0;
  //_beta = _r.norm(); pas besoin si beta existe déjà dans le main
}
void GMRes::Arnoldi(SparseVector<double> v, int N, SparseMatrix<double> A, SparseMatrix<double> v_arno, SparseMatrix<double>  H)
{
  //Définition des tailles des matrices H et v_arno, qui deviendront Hm et Vm
  v_arno.resize(N,N);
  H.resize(N,N);
  //definition des vecteurs et matrices locaux utilisés:
  SparseVector<double> wj(N); // Produit matrice vecteur A*vj, économiser des opérations
  SparseVector<double> vi(N); // Vecteur vi, économiser des opérations
  SparseVector<double> zj(N);
  SparseVector<double> Sk(N); // Vecteur résultant de la somme
  //Normalisation de v1
  v_arno.col(1)= (1/v.norm())*v;
  for (int j=1 ; j < N  ; j++)
  {
    wj= A*v_arno.col(j);
    for (int i=1 ; i < j  ; i++)
    {
      vi=v_ano.col(i);
      H(i,j)= wj.dot(vi);
    }
    for (int k=1 ; k < N  ; k++)
    {
      Sk = Sk + H(k,j)*v_arno.col(j);
    }
    zj = wj - Sk;
    H(j+1,j)= zj.norm();
    if (H(j+1,j)=0)
    {
      break;
    }
    v_arno.col(j+1) = zj / H(j+1,j) ;
  }
}
void GMRes::calcul_sol(SparseVector<double> b, SparseMatrix<double> A)
{}
*/

#define _METHODE_RES_CPP
#endif
