#include "ReadMatrix.h"
#include <iostream>
#include <fstream>


void InitialisationMatrixA(int& N,string name_Matrix, SparseMatrix<double>& A)
{
  double taille, nb_non_nul;
  //ouvrir le fichier nam_Matrix
  ifstream mon_flux(name_Matrix+".mtx");
  // recuperer la taille de la matrice et le nombre d'éléments non nuls
  string er1, er2, er3,er4, er5;
  mon_flux >> er1 >> er2 >> er3>> er4>> er5;
  //cout << er1 << ","<< er2 << "," << er3<<"," << er4<<"," << er5 << endl;
  mon_flux >> taille >> taille >> nb_non_nul;

  N=int(taille);
  //cout << N << ","<< taille << "," << nb_non_nul<< endl;
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
    //cout << l <<","<< valeur << endl;
    triplets.push_back({l-1,c-1,valeur});
if (er5== "symmetric" && l!=c)
{
  triplets.push_back({c-1,l-1,valeur});
}
  }
  mon_flux.close();
  A.setFromTriplets(triplets.begin(),triplets.end());
}
