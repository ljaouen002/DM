#ifndef _METHODE_RES_CPP

#include "MethodeRes.h"
#include <iostream>
#include "Dense"
#include "Sparse"
//#include<Eigen/SparseQR>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
MethodeRes::MethodeRes()
{}


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



    Jacobi::Jacobi(SparseMatrix<double> M, SparseMatrix<double> N)
    {
      _M=M;
      _N=N;
    }

    void Jacobi :: Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);


      //Définition des matrices D, E et F
      SparseMatrix<double> _D, _F, _E;
      _D=0*_A;
      _E=0*_A;
      _F=0*_A;
      _M=0*_A;
      _N=0*_A;

      //Création de M=D^-1
      for (int i=0 ; i<_D.rows() ; ++i)
      {
        _D.coeffRef(i,i) = _A.coeffRef(i,i);
        _M.coeffRef(i,i)=1/_D.coeffRef(i,i);
      }

      //Création de N=E+F
      _F = -_A.triangularView<StrictlyUpper>();
      _E = -_A.triangularView<StrictlyLower>();
      _N=_E+_F;

      _r=_b-_A*_sol0;

    }

    void Jacobi::calcul_sol(Eigen::SparseVector<double>& _r)
    {
      _sol=_M*_N*_sol+_M*_b;
      _r=_b-_A*_sol;
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
      //Définition des variables locales
      SparseVector<double> z;
      double _alpha;

      z = _A *_r;
      _alpha = _r.dot(_r)/z.dot(_r);
      _sol= _sol + _alpha*_r;
      _r=_r-_alpha*z;

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
    }






    Residu_Precondi_gauche:: Residu_Precondi_gauche(Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> F, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> D1, Eigen::SparseVector<double> q , int precondi)
    {
      _M=M;
      _precondi = precondi;
      _q=q;
      _F=F;
      _E=E;
      _D=D;
      _D1=D1;
    }

    void Residu_Precondi_gauche::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      _M.resize(_A.rows(), _A.cols());
      _q.resize(_A.rows());
      _F.resize(_A.rows(), _A.cols());
      _D.resize(_A.rows(), _A.cols());
      _E.resize(_A.rows(), _A.cols());
      _D1.resize(_A.rows(), _A.cols());

      _r=_b-_A*_sol0;


      if (_precondi == 1) //Jacobi
      {
        //Définition de M-1
        for (int i=0 ; i<_A.rows() ; ++i)
        {
          _M.coeffRef(i,i)=1/_A.coeffRef(i,i);
        }

        //Résolution
        _q= _M*_r;

      }




      else if (_precondi == 2) //SGS
      {

        double diago;
        SparseMatrix<double> Md;
        SparseVector<double> q1;
        Md.resize(_A.rows(), _A.cols());
        q1.resize(_A.rows());

        //Création de D, E et F
        for (int i=0 ; i< (_A.rows()) ; ++i)
        {
          _D.coeffRef(i,i) = _A.coeffRef(i,i);
          _D1.coeffRef(i,i)= 1/ _A.coeffRef(i,i);
        }

        _F = -_A.triangularView<StrictlyUpper>();
        _E = -_A.triangularView<StrictlyLower>();



        Md = (_D-_E)*_D1;

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = _r.coeffRef(i) - Md.row(i).dot(q1);
          q1.insert(i) = diago/ Md.coeffRef(i,i);
        }


        Md = (_D-_F);

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = q1.coeffRef(_A.rows() -1 -i) - Md.row(_A.rows() -1 -i).dot(_q);
          _q.insert(_A.rows() -1 -i) = diago/ Md.coeffRef(_A.rows() -1 -i, _A.rows() -1 -i);
        }




      }
      else
      {
        cout << "Ce préconditionnement n'existe pas!" <<endl ;
      }


    }

    void Residu_Precondi_gauche::calcul_sol(Eigen::SparseVector<double>& _r)
    {
      double _alpha;
      SparseVector<double> z, w;

      w=_A*_q;

      if (_precondi == 1) //Jacobi
      {

        //Résolution
        z=_M*w;
      }


      else if (_precondi == 2)

      {
        SparseMatrix<double> Md;
        double diago;
        SparseVector<double> z1;

        z1.resize(_A.rows());
        z.resize(_A.rows());


        Md = (_D-_E)*_D1;

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = w.coeffRef(i) - Md.row(i).dot(z1);
          z1.insert(i) = diago/ Md.coeffRef(i,i);
        }


        Md = (_D-_F);

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = z1.coeffRef(_A.rows() -1 -i) - Md.row(_A.rows() -1 -i).dot(z);
          z.insert(_A.rows() -1 -i) = diago/ Md.coeffRef(_A.rows() -1 -i, _A.rows() -1 -i);
        }

      }


      _alpha=(_q.dot(z))/(z.dot(z));
      cout << _alpha <<endl;
      _sol=_sol+_alpha*_q;
      _r=_r-_alpha*w;
      _q=_q- _alpha*z;
    }


    Residu_Precondi_droite:: Residu_Precondi_droite(Eigen::SparseMatrix<double> M, Eigen::SparseMatrix<double> F, Eigen::SparseMatrix<double> E, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> D1,  int precondi)
    {
      _M=M;
      _precondi = precondi;
      _F=F;
      _E=E;
      _D=D;
      _D1=D1;
    }

    void Residu_Precondi_droite::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      _M.resize(_A.rows(), _A.cols());
      _r=_b-_A*_sol0;

      if (_precondi == 3) //Jacobi
      {
        //Définition de M-1
        for (int i=0 ; i<_A.rows() ; ++i)
        {
          _M.coeffRef(i,i)=1/_A.coeffRef(i,i);
        }
      }

      else if (_precondi == 4)
      {

        _F.resize(_A.rows(), _A.cols());
        _D.resize(_A.rows(), _A.cols());
        _E.resize(_A.rows(), _A.cols());
        _D1.resize(_A.rows(), _A.cols());


        //Création de D, E et F
        for (int i=0 ; i< (_A.rows()) ; ++i)
        {
          _D.coeffRef(i,i) = _A.coeffRef(i,i);
          _D1.coeffRef(i,i)= 1/ _A.coeffRef(i,i);
        }

        _F = -_A.triangularView<StrictlyUpper>();
        _E = -_A.triangularView<StrictlyLower>();
      }

      else
      {
        cout << "Ce préconditionnement n'existe pas!" <<endl ;
      }

    }

    void Residu_Precondi_droite::calcul_sol(Eigen::SparseVector<double>& _r)
    {
      double _alpha;
      SparseVector<double> z, w;


      if (_precondi == 3)
      {
        z=_M*_r;
      }



      else if (_precondi == 4)
      {
        SparseMatrix<double> Md;
        double diago;
        SparseVector<double> z1;

        z1.resize(_A.rows());
        z.resize(_A.rows());


        Md = (_D-_E)*_D1;

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = _r.coeffRef(i) - Md.row(i).dot(z1);
          z1.insert(i) = diago/ Md.coeffRef(i,i);
        }


        Md = (_D-_F);

        for (int i = 0; i < _A.rows(); i++)
        {
          diago = z1.coeffRef(_A.rows() -1 -i) - Md.row(_A.rows() -1 -i).dot(z);
          z.insert(_A.rows() -1 -i) = diago/ Md.coeffRef(_A.rows() -1 -i, _A.rows() -1 -i);
        }
      }




      w=_A*z;
      _alpha=(_r.dot(w))/(w.dot(w));
      _sol=_sol+_alpha*z;
      _r=_r-_alpha*w;

    }




    Residu_Precondi_auto:: Residu_Precondi_auto()
    {

    }

    void Residu_Precondi_auto::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      _r=_b-_A*_sol0;

    }

    void Residu_Precondi_auto::calcul_sol(Eigen::SparseVector<double>& _r)

    {

      double _alpha;
      SparseVector<double> z, w;

      for (int i=0 ; i<300 ; ++i)
      {
        z=_A*_r;
        _alpha=(_r.dot(z))/(z.dot(z));
        _sol=_sol+_alpha*_r;
        _r=_r-_alpha*z;
      }


      w=_A*z;
      _alpha=(_r.dot(w))/(w.dot(w));
      _sol=_sol+_alpha*z;
      _r=_r-_alpha*w;
    }


    //Méthode de GMRes
    GMRes::GMRes(SparseMatrix<double> Vm, SparseMatrix<double> Hm, int m)
    {
      _Vm=Vm;
      _Hm=Hm;
      _m=m;
    }


    void GMRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);
      _r=_b-_A*_sol0;


    }
    void GMRes::Arnoldi(SparseVector<double> v, SparseMatrix<double> A, SparseMatrix<double> Vm, SparseMatrix<double> Hm)
    {

      SparseVector<double> wj(A.rows()); // Produit matrice vecteur A*vj, économiser des opérations
      SparseVector<double> vi(A.rows()); // Vecteur vi, économiser des opérations
      SparseVector<double> zj(A.rows());
      SparseVector<double> Sk(A.rows()); // Vecteur résultant de la somme

      _Vm.col(0)= v*(1/v.norm());

      for (int j=0 ; j < _m ; j++)
      {
        Sk.setZero();
        wj= A*_Vm.col(j);

        for (int i=0 ; i < j+1 ; i++)
        {
          vi=_Vm.col(i);
          _Hm.coeffRef(i,j)= wj.dot(vi);


          for (int k=0 ; k < _m ; k++)
          {
            Sk = Sk + _Hm.coeffRef(k,j)*_Vm.col(j);
          }

          zj = wj - Sk;
        }


        _Hm.coeffRef(j+1,j)= zj.norm();


        if (_Hm.coeff(j+1,j)< 1e-14)
        {
          _Hm.coeffRef(j+1,j)=0;
          cout << "IF LE COEFF  Hj+1j  EST NUL  " <<  endl;
          break;
        }
        _Vm.col(j+1) = zj / _Hm.coeffRef(j+1,j) ;

      }
      //Verification de l'orthogonalité des colonnes de la matrice Vm
      //  cout << _Vm.col(1).dot(_Vm.col(0)) << endl;


    }




    void GMRes::calcul_sol(Eigen::SparseVector<double>& _r)
    {

      SparseVector<double> e1, gm, y2, y;
      double diago;
            double beta;


      //Premier vecteur de la base canonique
      e1.resize(_m+1);
      e1=0*e1;
      e1.coeffRef(0)=1;

      //Solution équation normale
      y.resize(_m);
      y2.resize(_m+1);

      //vecteur gm
      gm.resize(_m+1);

      //Matrice obtenue par décomprisaition QR
      SparseMatrix<double> Qm, Rm;
      _Hm.resize(_m+1,_m);
      _Vm.resize(_A.rows(),_m+1);
      Qm.resize(_m+1,_m+1);
      Rm.resize(_m+1,_m);

            SparseMatrix<double> Vmy;


      //Application d'arnoldi à r
      GMRes::Arnoldi(_r, _A, _Vm, _Hm);


      beta=_r.norm();

      //Décompostion QR
      Eigen::SparseQR< SparseMatrix<double>, COLAMDOrdering<int> > solver_direct;
      _Hm.makeCompressed();
      solver_direct.compute(_Hm);
      Qm=solver_direct.matrixQ();
      Rm=solver_direct.matrixR();


      gm= beta*Qm.transpose()*e1;


      y=0*y;

      for (int i = 0; i < _m; i++)
      {
        diago = gm.coeffRef(_m-1-i) - Rm.row(_m-1-i).dot(y);
        y.insert(_m-1-i) = diago/Rm.coeffRef(_m-1-i,_m-1-i);
      }



      for (int i = 0; i < _m; i++)
      {
        y2.coeffRef(i)=y.coeffRef(i);
      }
      y2.coeffRef(_m)=0;


      y.resize(_m+1);
      y=y2;

          Vmy = _Vm*y;

      _sol= _sol + Vmy;
      _r= _r - _A*Vmy;
      beta=_r.norm();

    }

    #define _METHODE_RES_CPP
    #endif
