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





    // void MethodeRes::InitialisationMat(int& N,string name_Matrix, SparseMatrix<double>& _A)
    // {
    //   double taille, nb_non_nul;
    //   //ouvrir le fichier nam_Matrix
    //   ifstream mon_flux(name_Matrix+".mtx");
    //   // recuperer la taille de la matrice et le nombre d'éléments non nuls
    //   string er1, er2, er3,er4, er5;
    //   mon_flux >> er1 >> er2 >> er3>> er4>> er5;
    //   //cout << er1 << ","<< er2 << "," << er3<<"," << er4<<"," << er5 << endl;
    //   mon_flux >> taille >> taille >> nb_non_nul;
    //
    //   N=int(taille);
    //   //cout << N << ","<< taille << "," << nb_non_nul<< endl;
    //   _A.resize(N,N);
    //   // indice et valeur pour former la matrice
    //   int l, c;
    //   double valeur;
    //   vector<Triplet<double>> triplets;
    //   for (int k=0; k<int(nb_non_nul) ; k++)
    //   {
    //     mon_flux >> l;
    //     mon_flux >> c;
    //     mon_flux >> valeur;
    //     //cout << l <<","<< valeur << endl;
    //     triplets.push_back({l-1,c-1,valeur});
    //   }
    //   mon_flux.close();
    //   _A.setFromTriplets(triplets.begin(),triplets.end());
    // }

    // Initialisation de la matrices
    // void MethodeRes::InitialisationMat(string name_Matrix, SparseMatrix<double> _A, SparseVector<double> _b, SparseVector<double> _sol0, SparseVector<double>& _r)
    // {
    //   double taille, nb_non_nul;
    //   //ouvrir le fichier nam_Matrix
    //   ifstream mon_flux(name_Matrix+".mtx");
    //   // recuperer la taille de la matrice et le nombre d'éléments non nuls
    //   string er1, er2, er3,er4, er5;
    //   mon_flux >> er1 >> er2 >> er3>> er4>> er5;
    //   //cout << er1 << ","<< er2 << "," << er3<<"," << er4<<"," << er5 << endl;
    //   mon_flux >> taille >> taille >> nb_non_nul;
    //   int N;
    //   N=int(taille);
    //   //cout << N << ","<< taille << "," << nb_non_nul<< endl;
    //   _A.resize(N,N);
    //   // indice et valeur pour former la matrice
    //   int l, c;
    //   double valeur;
    //   vector<Triplet<double>> triplets;
    //   for (int k=0; k<int(nb_non_nul) ; k++)
    //   {
    //     mon_flux >> l;
    //     mon_flux >> c;
    //     mon_flux >> valeur;
    //     //cout << l <<","<< valeur << endl;
    //     triplets.push_back({l-1,c-1,valeur});
    //   }
    //   mon_flux.close();
    //
    //   _A.setFromTriplets(triplets.begin(),triplets.end());
    //
    //
    //   // Définition des vecteurs sol0 et b
    //   _sol0.resize(N);
    //   _b.resize(N);
    //   cout << "ici"<< endl;
    //   for (int i=0 ; i<_sol0.rows() ; i++)
    //   {
    //     _sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
    //     _b.coeffRef(i)=1.;          //Définir un valeur de b
    //   }
    // }


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
      //Résolution de Jacobi
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






    Residu_Precondi_gauche:: Residu_Precondi_gauche(Eigen::SparseMatrix<double> M, Eigen::SparseVector<double> q , int precondi)
    {
      _M=M;
      _precondi = precondi;
      _q=q;
    }

    void Residu_Precondi_gauche::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      _M.resize(_A.rows(), _A.cols());
      _q.resize(_A.rows());
      _r=_b-_A*_sol0;

      if (_precondi == 1) //Jacobi
      {
        //Définition de M-1
        for (int i=0 ; i<_A.rows() ; ++i)
        {
          _M.coeffRef(i,i)=1/_A.coeffRef(i,i);
        }

      }




      else if (_precondi == 2) //SGS
      {
        MatrixXd Md;
        Md.resize(_A.rows(), _A.cols());

        MatrixXd D, E, F, D1;
        F.resize(_A.rows(), _A.cols());
        D.resize(_A.rows(), _A.cols());
        E.resize(_A.rows(), _A.cols());
        D1.resize(_A.rows(), _A.cols());

        //Création de D, E et F
        for (int i=0 ; i< (_A.rows()) ; ++i)
        {
          D(i,i) = _A.coeffRef(i,i);
          D1(i,i)= 1/ _A.coeffRef(i,i);
        }

        F = -_A.triangularView<StrictlyUpper>();
        E = -_A.triangularView<StrictlyLower>();

        //Calcul de M
        Md = (D-E)*D1*(D-F);

        //Inversion de M
        Md= Md.inverse();
        _M=Md.sparseView();
      }
      else
      {
        cout << "Ce préconditionnement n'existe pas!" <<endl ;
      }

      //Résolution
      _q= _M*_r;

    }

    void Residu_Precondi_gauche::calcul_sol(Eigen::SparseVector<double>& _r)
    {
      double _alpha;
      SparseVector<double> z, w;

      w=_A*_q;

      //Résolution
      z=_M*w;

      _alpha=(_q.dot(z))/(z.dot(z));
      _sol=_sol+_alpha*_q;
      _r=_r-_alpha*w;
      _q=_q- _alpha*z;
    }


    Residu_Precondi_droite:: Residu_Precondi_droite(Eigen::SparseMatrix<double> M,  int precondi)
    {
      _M=M;
      _precondi = precondi;
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


      else if (_precondi == 4) //SGS
      {
        MatrixXd Md;
        Md.resize(_A.rows(), _A.cols());

        MatrixXd D, E, F, D1;
        F.resize(_A.rows(), _A.cols());
        D.resize(_A.rows(), _A.cols());
        E.resize(_A.rows(), _A.cols());
        D1.resize(_A.rows(), _A.cols());

        //Création de D, E et F
        for (int i=0 ; i< (_A.rows()) ; ++i)
        {
          D(i,i) = _A.coeffRef(i,i);
          D1(i,i)= 1/ _A.coeffRef(i,i);
        }

        F = -_A.triangularView<StrictlyUpper>();
        E = -_A.triangularView<StrictlyLower>();

        //Calcul de M
        Md = (D-E)*D1*(D-F);

        Md= Md.inverse();
        _M=Md.sparseView();
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

      //Résolution
      z=_M*_r;
      w=_A*z;
      _alpha=(_r.dot(w))/(w.dot(w));
      _sol=_sol+_alpha*z;
      _r=_r-_alpha*w;

    }




    Residu_Precondi_auto:: Residu_Precondi_auto(Eigen::SparseMatrix<double> M)
    {
      _M=M;
    }

    void Residu_Precondi_auto::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
      MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      // MatrixXd Md;
      // Md.resize(_A.rows(), _A.cols());
      // _M.resize(_A.rows(), _A.cols());
      _r=_b-_A*_sol0;
      //
      // for (int i=0 ; i<A.rows() ; ++i)
      // {
      //   for (int j=0 ; j<A.cols() ; ++j)
      //   {
      //     Md(i,j)=A.coeffRef(i,j);
      //   }
      // }
      //
      // Md= Md.inverse();
      // _M=Md.sparseView();


    }

    void Residu_Precondi_auto::calcul_sol(Eigen::SparseVector<double>& _r)

    {

      double _alpha;
      SparseVector<double> z, w;
      //Résolution

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

        }


        for (int k=0 ; k < _m ; k++)
        {
          Sk = Sk + _Hm.coeffRef(k,j)*_Vm.col(j);
        }


        zj = wj - Sk;


        _Hm.coeffRef(j+1,j)= zj.norm();


        if (_Hm.coeff(j+1,j)< 1e-14)
        {
          _Hm.coeffRef(j+1,j)=0;
          cout << "IF LE COEFF  Hj+1j  EST NUL  " <<  endl;
          break;
        }
        _Vm.col(j+1) = zj / _Hm.coeffRef(j+1,j) ;

      }

      //  cout << _Vm.col(1).dot(_Vm.col(0)) << endl;

      //Eigen::VectorXd _w;
      // _Hm.resize(_m+1,_m);
      // _Vm.resize(_A.rows(),_m+1);
      // //v=_r;
      // Eigen :: SparseVector<double> _w;
      // _w.resize(_A.cols());
      // _w.setZero();
      // _Vm.col(0) = v/v.norm();
      //
      // //--------------------------------------------------------------------
      // //------------------------BOUCLE SUR K------------------------
      // //------------------------------------------------------------------------
      // for (int k=0; k< _m; k++)
      // {
      //   _w = _A*_Vm.col(k);
      //
      // //_____________________________________________________
      // //___________________BOUCLE SUR I_____________________
      // //____________________________________________________
      //   for (int i=0; i < k+1; i++)
      //   {
      //     _Hm.coeffRef(i,k)=_w.dot(_Vm.col(i));
      //     _w = _w-_Hm.coeffRef(i,k)*_Vm.col(i);
      //   }
      // //----------------------------------------------------------------------------
      // //-------------------ON SORT DE LA BOUCLE EN I--------------------------------
      // //-----------------------------------------------------------------------------
      //
      //   _Hm.coeffRef(k+1,k) = _w.norm();
      //
      //   if (_Hm.coeffRef(k+1,k) ==0)
      //   {
      //     break;
      //   }
      //
      //   _Vm.col(k+1)=_w/_Hm.coeffRef(k+1,k);
      //
      // }
      // cout << "_____________Vm__________________ " <<endl;
      // cout << _Vm << endl;
      // cout <<"__________hmmmmmmm_____________"<<endl;
      // cout << _Hm<<endl;
      // cout << "-------------V0.V1-----------" << endl;
      // cout << _Vm.col(0).dot(_Vm.col(1)) << endl;


    }




    void GMRes::calcul_sol(Eigen::SparseVector<double>& _r)
    {

      SparseVector<double> e1, gm, y2;
      double b;


      VectorXd y;

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
      SparseMatrix<double> Qm, Rm, AVm;

      _Hm.resize(_m+1,_m);
      _Vm.resize(_A.rows(),_m+1);
         AVm.resize(_A.rows(),_m+1);
      Qm.resize(_m+1,_m+1);
      Rm.resize(_m+1,_m);

      AVm= _A*_Vm;

      //Application d'arnoldi à r
      GMRes::Arnoldi(_r, _A, _Vm, _Hm);
      //  cout << "Hm" << _Hm << endl;

      double beta;
      beta=_r.norm();

      //Décompostion QR
      //template<typename _MatrixType , typename _OrderingType >
      Eigen::SparseQR< SparseMatrix<double>, COLAMDOrdering<int> > solver_direct;
      //Pour utiliser cette fonction, doit compresser Hm
      _Hm.makeCompressed();
      solver_direct.compute(_Hm);
      Qm=solver_direct.matrixQ();
      Rm=solver_direct.matrixR();

      //
      // cout << "Hm" << _Hm << endl;
      // cout << "prod" << Qm*Rm << endl;

      //    gm= beta*Qm.transpose()*e1;
      gm= beta*e1;

      //
      y=0*y;
      //  y.coeffRef(0)=2.65409

    //  solver_direct.compute(AVm);
    //


    //  y=  solver_direct.solve(_r);


      //  Qte1 = _Beta*Qte1;


      //
      // for (int i = 0; i < _m; i++)
      // {
      //   b = gm.coeffRef(_m-1-i) -_Hm.row(_m-1-i).dot(y);
      //   //cout << i << "H      " << _Hm.row(_m-1-i)  << endl;
      //   //cout << i << "g      " << gm.coeffRef(_m-1-i)  << endl;
      // //  cout << i << "b      " << b  << endl;
      //   y.insert(_m-1-i) = b/_Hm.coeffRef(_m-1-i,_m-1-i);
      // //  cout << i << y  << endl;
      // }
      //
      //
      //
      // for (int i = 0; i < _m; i++)
      // {
      //   y2.coeffRef(i)=y.coeffRef(i);
      // }
      // y2.coeffRef(_m)=0;
      //
      //
      // y.resize(_m+1);
      // y=y2;
      //
      //
      // _sol= _sol + _Vm*y;
      // _r= _r-_A*_Vm*y;
      // beta=_r.norm();

    }

    #define _METHODE_RES_CPP
    #endif

