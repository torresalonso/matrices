/*Programa que lee un polinomio y encuentra una raíz utilizando el método de Bisección.
Este programa utiliza sólo tipo de dato long double para mostrar su funcionamiento en casos extremos.*/
#include<iostream>
#include "matrices.h"
#include<math.h>

using namespace std;

void ClearScreen();

int main(){

  int n=0, i=0,j=0;
  char normalizar = 'S';

  cout << "n:" << '\n';
  cin>>n;

  //matrices cuadradas
  matriz_cuadrada A, aux;
  A.dimension = n;
  aux.dimension = n;

  //leer la matriz
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
      cout<<"A_"<<i+1<<","<<j+1<<"\t";
      cin >> A.matriz[i][j];
    }

  ClearScreen();

  cout<<"\n";
  //-Imprimir A
  cout << "--A:\n";
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout << A.matriz[i][j] << "\t";
    cout<<"\n";
  }

  cout<<"\n";
  //-transpuesta
  cout<<"--transpuesta\n";
  aux = matriz_transpuesta(A);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout<<aux.matriz[i][j]<<"\t";
    cout<<"\n";
  }

  cout<<"\n";
  //-producto
  cout<<"--producto A²\n";
  aux = producto_matrices(A, A);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout<<aux.matriz[i][j]<<"\t";
    cout<<"\n";
  }

  cout<<"\n";
  //-ortogonal Q
  cout<<"--Ortogonalizar las columnas con gram_schmidt\n";
  matriz_cuadrada Q;
  Q.dimension = A.dimension;
  cout<<"ortonormal? S/n: \n";
  cin >> normalizar;
  if(normalizar == 's' || normalizar == 'S')
    Q = QR_getQ(A, 1);
  else
    Q = QR_getQ(A, 0);

ClearScreen();

  cout<<"--Q:\n";
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout<<Q.matriz[i][j]<<"\t";
    cout<<"\n";
  }
  cout<<"\n";

  matriz_cuadrada R;
  R.dimension = A.dimension;
  R = QR_getR(A, Q);

  cout<<"--R:\n";
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout<<R.matriz[i][j]<<"\t";
    cout<<"\n";
  }

  cout<<"\n";
  //-Algoritmo QR
  cout<<"--Algoritmo QR\n";
  A = QR_eigenvalores(A, 50);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout << A.matriz[i][j]<<"\t";
    cout<<"\n";
  }
  return 0;
}

void ClearScreen()
{ cout << string( 100, '\n' );}
