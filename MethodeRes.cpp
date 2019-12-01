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
//SparseMatrix<double> D, SparseMatrix<double> F, SparseMatrix<double> E,

Jacobi::Jacobi(SparseMatrix<double> M, SparseMatrix<double> N)
{
  //_D=D;
//  _F=F;
  //_E=E;
  _M=M;
  _N=N;
}

void Jacobi :: Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)

{
  MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

SparseMatrix<double> _D, _F, _E;

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
  // for (int i=0 ; i<_A.rows() ; ++i)
  // {
  //   for (int j=0 ; j<_A.cols() ; ++j)
  //   {
  //     if (i<j)
  //     {
  //       _E.coeffRef(i,j)=-_A.coeffRef(i,j);     // Partie triangulaire supérieure de A
  //     //  cout << "E"<<_E.coeffRef(i,j) << endl;
  //     }
  //     else if (i>j)
  //     {
  //       _F.coeffRef(i,j)=-_A.coeffRef(i,j);     // Partie triangulaire inférieure de A
  //     //  cout << "F" << _F.coeffRef(i,j) << endl;
  //     }
  //   }
  // }

  // F= - triangle inférieur de A
_F = - _A.triangularView<StrictlyUpper>();
// F= - triangle supérieur de A
_E = - _A.triangularView<StrictlyLower>();
  _N=_E+_F;

  _r=_b-_A*_sol0;
 /*cout << "M" << _M  << endl;
cout << "N" << _N << endl;
cout << "_r=" << _r << endl;
  cout << "rcalc" << _r.norm() << endl;*/

}


void Jacobi::calcul_sol(Eigen::SparseVector<double>& _r)
{


  _sol=_M*_N*_sol+_M*_b;
  cout << "Sol" << _sol << endl;

  _r=_b-_A*_sol;
//  cout <<_N.norm() << " " <<  "1" << endl;
//  cout <<_r.norm() << " " <<  "2" << endl;
//cout << "rcalc" << _r.norm() << endl;


}



//Méthode du GPO

GPO::GPO()
{

}


void GPO :: Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
    MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

    _r=_b-_A*_sol0;
      cout << "r" << _r << endl;
}

void GPO :: calcul_sol(Eigen::SparseVector<double>& _r)
{
  SparseVector<double> z;
  double _alpha;

  z = _A *_r;
  _alpha = _r.dot(_r)/z.dot(_r);
  _sol= _sol + _alpha*_r;
 _r=_r-_alpha*z;


 cout << "z" << z  << endl;
 cout << "alpha" << _alpha << endl;
 cout << "_sol" << _sol << endl;
  cout << "r" << _r << endl;
}


//Méthode du résidu minimum

Residu:: Residu()
{

}

void Residu::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
  MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

  _r=_b-_A*_sol0;
}

void Residu::calcul_sol(Eigen::SparseVector<double>& _r)
{
  double _alpha;
  SparseVector<double> _z;



  _z=_A*_r;
  _alpha=_r.dot(_z)/_z.dot(_z);
  _sol=_sol+_alpha*_r;
  _r=_r-_alpha*_z;
    cout << "r" << _r << endl;
}




//Méthode de GMRes
GMRes::GMRes(SparseMatrix<double> v_arno, SparseMatrix<double> H)
{
  _v_arno=v_arno;
  _H=H;
}
void GMRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
{
  _r=b-A*_sol0;

}
void GMRes::Arnoldi(SparseVector<double> v, SparseMatrix<double> A, SparseMatrix<double> v_arno, SparseMatrix<double>  H)
{
  // //Définition des tailles des matrices H et v_arno, qui deviendront Hm et Vm
  // //_v_arno.resize(A.rows(),A.cols());
  //_H.resize(A.rows(),A.cols());
  // //definition des vecteurs et matrices locaux utilisés:
  SparseVector<double> wj(A.rows()); // Produit matrice vecteur A*vj, économiser des opérations
  SparseVector<double> vi(A.rows()); // Vecteur vi, économiser des opérations
  SparseVector<double> zj(A.rows());
  SparseVector<double> Sk(A.rows()); // Vecteur résultant de la somme
//  MatrixXd H1(A.rows(),A.cols())
//  double g;
  //Normalisation de v1
 v_arno.col(1)= (1/v.norm())*v;
  for (int j=1 ; j < A.rows() ; j++)
 {
    wj= A*v_arno.col(j);
    for (int i=1 ; i < j  ; i++)
     {
      vi=v_arno.col(i);
       H.coeffRef(i,j)= wj.dot(vi);
    }
    for (int k=1 ; k < A.rows()  ; k++)
    {
      Sk = Sk + H.coeffRef(k,j)*v_arno.col(j);
    }
    zj = wj - Sk;
    H.coeffRef(j+1,j)= zj.norm();

    if (H.coeffRef(j+1,j)=0)
    {
      break;
   }
    v_arno.col(j+1) = zj / H.coeffRef(j+1,j) ;
}
}



void GMRes::calcul_sol(Eigen::SparseVector<double>& _r)
{
  SparseMatrix<double> Hm, Vm;
  GMRes::Arnoldi(_r, _A, Vm, Hm);

  
}

#define _METHODE_RES_CPP
#endif
