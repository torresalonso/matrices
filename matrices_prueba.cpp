/*Programa que lee un polinomio y encuentra una raíz utilizando el método de Bisección.
Este programa utiliza sólo tipo de dato long double para mostrar su funcionamiento en casos extremos.*/
#include<iostream>
#include "matrices.h"

using namespace std;

void ClearScreen();

int main(){

  int n=0, i=0,j=0;
  char normalizar = 'S';

  cout << "n:" << '\n';
  cin>>n;
  cout<<"ortonormal? S/n: \n";
  cin>>normalizar;

  //matrices cuadradas
  matriz_cuadrada A;
  A.dimension=n;

  //leer la matriz
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
      cout<<"A_"<<i+1<<","<<j+1<<"\t";
      cin >> A.matriz[i][j];
    }

ClearScreen();

  matriz_cuadrada Q;
  Q=gram_schmidt(A);

  //imprimir matriz ortogonalizada
  cout<<"\nMatriz ortonormal\n";
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout << Q.matriz[i][j] << "\t";
    cout << "\n";
  }

  return 0;
}

void ClearScreen()
{ cout << string( 100, '\n' );}
