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



    // Initialisation de la matrices
    void MethodeRes::InitialisationMat(string name_Matrix, SparseMatrix<double> A, SparseVector<double> b, SparseVector<double> sol0, SparseVector<double> r)
    {
      double taille, nb_non_nul;
      //ouvrir le fichier nam_Matrix
      ifstream mon_flux(name_Matrix+".mtx");
      // recuperer la taille de la matrice et le nombre d'éléments non nuls
      mon_flux >> taille >> taille >> nb_non_nul;
      int N;
      N=int(taille);
      A.resize(N,N);
      // indice et valeur pour former la matrice
      int l, c;
      double valeur;
      vector<Triplet<double>> triplets;
      for (int k=0; k<int(nb_non_nul) ; k++)
      {
        mon_flux >> l;
        mon_flux >> c;
        mon_flux >> valeur;
        triplets.push_back({l-1,c-1,valeur});
      }
      mon_flux.close();
      A.setFromTriplets(triplets.begin(),triplets.end());

      // Définition des vecteurs sol0 et b
    	sol0.resize(N) ; b.resize(N) ; r.resize(N);
    	for (int i=0 ; i<sol0.rows() ; i++)
    	{
    		sol0.coeffRef(i)=1.;       //Définir un valeur de sol0
    		b.coeffRef(i)=1.;          //Définir un valeur de b
    	}
    }

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
      //cout << "r" << _r << endl;
    }



    //Méthode de GMRes
    GMRes::GMRes(SparseMatrix<double> v_arno, SparseMatrix<double> H)
    {
      _v_arno=v_arno;
      _H=H;
    }
    void GMRes::Initialisation(SparseVector<double> b, SparseMatrix<double> A, SparseVector<double> sol0, Eigen::SparseVector<double>& _r,string results, MethodeRes* methode)
    {
            MethodeRes::Initialisation(b,A,sol0,_r,results,methode);

      _r=_b-_A*_sol0;

    }
    void GMRes::Arnoldi(SparseVector<double> v, SparseMatrix<double> A, SparseMatrix<double> v_arno, SparseMatrix<double> H)
    {
      // //Définition des tailles des matrices H et v_arno, qui deviendront Hm et Vm
      //_v_arno.resize(_A.rows(),_A.cols());
      //_H.resize(A.rows(),A.cols());
      // //definition des vecteurs et matrices locaux utilisés:
      SparseVector<double> wj(A.rows()); // Produit matrice vecteur A*vj, économiser des opérations
      SparseVector<double> vi(A.rows()); // Vecteur vi, économiser des opérations
      SparseVector<double> zj(A.rows());
      SparseVector<double> Sk(A.rows()); // Vecteur résultant de la somme
      // MatrixXd H1(A.rows(),A.cols())
      // double g;
      //Normalisation de v1
      v_arno.col(1)= (1/v.norm())*v;
      for (int j=1 ; j < A.rows()-1 ; j++)
      {
        wj= A*v_arno.col(j);
        for (int i=1 ; i < j+1 ; i++)
        {
          vi=v_arno.col(i);
          H.coeffRef(i,j)= wj.dot(vi);
        }
        for (int k=1 ; k < A.rows() ; k++)
        {
          Sk = Sk + H.coeffRef(k,j)*v_arno.col(j);
        }
        zj = wj - Sk;
        H.coeffRef(j+1,j)= zj.norm();

        if (H.coeffRef(j+1,j)==0)
        {
          break;
        }
        v_arno.col(j+1) = zj / H.coeffRef(j+1,j) ;
      }
    }





    void GMRes::calcul_sol(Eigen::SparseVector<double>& _r)
    {

      SparseVector<double> e1, y, gm;

      //Premier vecteur de la base canonique
      e1.resize(_A.rows()+1);
      e1=0*e1;
      e1.coeffRef(0)=1;

      //Solution équation normale associée à AVmy
      y.resize(_A.rows()+1);

      //vecteur gm
      gm.resize(_A.rows()+1);

      //Matrice obtenue par Arnoldi
      SparseMatrix<double> Hm, Vm;

      //Matrice obtenue par décomprisaition QR
      SparseMatrix<double> Qm, Rm;

      Hm.resize(_A.rows()+1,_A.cols());
      Vm.resize(_A.rows(),_A.cols());

      Qm.resize(_A.rows()+1,_A.cols()+1);
      Rm.resize(_A.rows()+1,_A.cols());



      //Application d'arnoldi à r
      GMRes::Arnoldi(_r, _A, Vm, Hm);

      double beta;
      beta=_r.norm();

      //Décompostion QR
      //template<typename _MatrixType , typename _OrderingType >
      Eigen::SparseQR< SparseMatrix<double>, COLAMDOrdering<int> > solver_direct;
      //Pour utiliser cette fonction, doit compresser Hm
      Hm.makeCompressed();
      solver_direct.compute(Hm);
      Qm=solver_direct.matrixQ();
      Rm=solver_direct.matrixR();

      gm= beta*Qm.transpose()*e1;
      //On souhaite que Rmy=gm, on retrouve la forme Ax=b
      Rm.makeCompressed();
      solver_direct.compute(Rm);
      y=solver_direct.solve(gm);


      _sol= _sol + Vm*y;
      _r= beta*e1 -Hm*y;
      beta=_r.norm();


      //     // Décompostion QR de Hm
      //     SparseQR<SparseMatrix<double>,COLAMDOrdering<int>> solver;
      //     Hm.makeCompressed();
      //     solver.compute(Hm);
      //     Qm=solver.matrixQ();
      //     Rm=solver.matrixR();
      //     gm=beta1*Qm.transpose()*e1;
      //     solver.compute(Rm);
      //     // Résolution de Rm*y=gm
      //     y=solver.solve(gm);
      //
      //     _X=_X+Vm*y;
      //     // calcul du nouveau résidu
      //     _r=_r-_A*Vm*y;


      // {
      //   // Création des différentes matrices et vecteurs utiles dans cette méthode
      //   SparseMatrix<double> Hm,Qm,Rm;
      //   SparseMatrix<double> Vm;
      //   VectorXd y,e1,gm;
      //   SparseVector<double> zj,a,r;
      //   double beta1;
      //   int k;
      //
      //   // Allocation des tailles des matrices et vecteurs
      //   r.resize(_N);
      //   zj.resize(_N);
      //   a.resize(_N);
      //   Vm.resize(_N,_m);
      //   Hm.resize(_m+1,_m);
      //   Qm.resize(_m+1,_m+1);
      //   gm.resize(_m+1);
      //   Rm.resize(_m+1,_m);
      //   y.resize(_m);
      //   e1.resize(_m+1);
      //
      //   // Création du m-ième vecteur de la base canonique de taille m
      //   for (int i=1; i<e1.size(); i++)
      //   {
      //     e1.coeffRef(i)=0.;
      //   }
      //   e1.coeffRef(0)=1.;
      //
      //   beta1=_r.norm();
      //   k=1;
      //   // Creation et ouverture du fichier GMRES.txt
      //   _file_out.open("GMRes"+to_string(_N)+".txt");
      //     _file_out<< k << "  " << beta1 << endl;
      //
      //   // Tant que la tolérence n'est pas atteinte ou que le nombre d'itérations maximum n'est pas réalisé
      //   while ((beta1>eps) && (k<kmax+1))
      //   {
      //     // Réalise une itération d'Arnoldi
      //     r=_r.sparseView();
      //     Vm.col(0)=r/r.norm();
      //     for (int j = 0; j < _m-1; j++)
      //     {
      //       a.setZero();
      //       for (int i = 0; i < j+1; i++)
      //       {
      //         Hm.coeffRef(i,j)=(_A*Vm.col(j)).dot(Vm.col(i));
      //         a+=Hm.coeff(i,j)*Vm.col(i);
      //       }
      //       zj=(_A*Vm.col(j))-a;
      //       Hm.coeffRef(j+1,j)=zj.norm();
      //       if (Hm.coeff(j+1,j)<pow(10,-14))
      //       {
      //         Hm.coeffRef(j+1,j)=0;
      //         break;}
      //       Vm.col(j+1)=zj/Hm.coeff(j+1,j);
      //     }
      //     // Fin de Arnoldi
      //
      //     // Décompostion QR de Hm
      //     SparseQR<SparseMatrix<double>,COLAMDOrdering<int>> solver;
      //     Hm.makeCompressed();
      //     solver.compute(Hm);
      //     Qm=solver.matrixQ();
      //     Rm=solver.matrixR();
      //     gm=beta1*Qm.transpose()*e1;
      //     solver.compute(Rm);
      //     // Résolution de Rm*y=gm
      //     y=solver.solve(gm);
      //
      //     _X=_X+Vm*y;
      //     // calcul du nouveau résidu
      //     _r=_r-_A*Vm*y;
      //
      //     beta1=_r.norm();
      //     k+=1;
      //
      //     // écriture dans le fichier
      //     _file_out<< k << "  " << beta1 << endl;
      //   }
      //   _file_out.close();
      //   if (k>kmax)
      //     {
      //       cout << "Tolérance non atteinte" << _r.norm() << endl;
      //     }
      // }


    }

    #define _METHODE_RES_CPP
    #endif
