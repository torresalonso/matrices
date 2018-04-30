/*Programa que lee un polinomio y encuentra una raíz utilizando el método de Bisección.
Este programa utiliza sólo tipo de dato long double para mostrar su funcionamiento en casos extremos.*/
#include<iostream>
#include<iomanip>
#include<math.h>
#include <unistd.h>

using namespace std;

double prodint(int longitud, double A[], double B[]);
double normalizar_vector(int longitud, double A[]);
void ClearScreen();
int main(){
  cout.precision(6);

  int n=0, i=0,j=0;
  char normalizar = 'S';

  cout << "n:" << '\n';
  cin>>n;
  cout<<"ortonormal? S/n: \n";
  cin>>normalizar;

  double a[n], q[n], aux[n];
//leer la matriz
  double A[n][n];
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
      cout<<"A_"<<i+1<<","<<j+1<<"\t";
      cin >> A[i][j];
    }

ClearScreen();

  //imprimir matriz
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      cout << A[i][j] << "\t";
    cout << "\n";
  }


  usleep((unsigned int) 50000);

  //-- calcular normas:
  for(j = 0; j < n; j++){
    for(i = 0; i < n; i++)
      a[i] = A[i][j];
    cout<<"||C"<<j+1<<"|| = " << sqrt(prodint(n,a,a))<<"\n";
  }

  int k=0, l=0, flag=0;

  //-- calcular ángulo entre cada par de columnas:
  for(k = 0; k+1 < n; k++){

    //copiar la columna k-ésima  en A
    for(i = 0; i < n; i++)
      a[i] = A[i][k];

    //este ciclo nos mueve por las columnas a la derecha de la columna k-ésima
    for(l = k+1; l < n; l++){

      //copiar en B la columna l-ésima (a la derecha de la columna k)
      for(i = 0; i < n ;i++)
        aux[i] = A[i][l];

      //enviar a calcular el producto interior
      cout << "<C" << k+1 << " , C" << l+1 << "> = " << prodint(n,a,aux)<<"\n";
      if(prodint(n,a,aux) != 0)
        flag++;
    }

  }

  if (flag == 0)
    cout << "Los vectores columna de la matriz son ORTO\n";
  else
    cout << "Hay " << flag << " vectores no ORTO\n";

  //---------------gram-schmidt---------------
  //iniciar la matriz ortogonal y el vector auxiliar
  double Q[n][n];
  for(i = 0; i < n; i++){
    aux[i] = 0;
    for (j = 0; j < n; j++)
      Q[i][j] = 0;
  }

  //q1 = a1
  for(i = 0; i < n; i++)
    Q[i][0] = q[i] = A[i][0];
  if(normalizar == 'S' || normalizar == 's')
    for(i = 0; i < n; i++)
      Q[i][0] = A[i][0]/sqrt(prodint(n,q,q));

  //--Obtener los vectores q_j, para j>1:
  /*el vector q almacena las columnas ya ortogonalizadas*/
  for(j = 1; j < n; j++){

    //tomar la j-ésima columna de la matriz A
    for(i = 0; i < n; i++)
      a[i] = A[i][j];

    //reinicia el vector auxiliar
    for(i = 0; i < n; i++)
      aux[i] = 0;

    //tomar las columnas desde la 0 hasta la (j-1)-ésima de la matriz ortogonal
    for(k = 0; k < j; k++){
      //copiar el vector q
      for(i = 0; i < n; i++)
        q[i]=Q[i][k];

      //el vector aux guardará la combinacion lineal
      for(i = 0; i < n; i++)
        aux[i] += ((prodint(n,a,q))/prodint(n,q,q))*q[i];
    }//al cabo de este for, ya tenemos la combinación lineal de los vectores q

    //Ortogonalizamos la j-ésima columna de A:
    for(i=0;i<n;i++)
      Q[i][j] = q[i] = a[i]-aux[i];

    //Normalizamos el vector q
    if(normalizar=='s'||normalizar=='S'){
      for(i=0;i<n;i++)
        Q[i][j] /= sqrt(prodint(n,q,q));
    }

      /*imprimir la matriz ortogonal*/
      for(i = 0; i < n; i++){
        for (int h = 0; h < n; h++)
          cout << Q[i][h] << "\t";
        cout << "\n\n";
      }


/*
      if(normalizar=='s'||normalizar=='S'){
        for(j=0;j<n;j++){
          for(i=0;i<n;i++)
            q[i] = (a[i]-aux[i]);
          for(i=0;i<n;i++)
            Q[i][j] /= prodint(n,q,q);
        }

        //imprimir la matriz ortonormal//
        cout<<"normalizada: \n";
        for(i = 0; i < n; i++){
          for (j = 0; j < n; j++)
            cout << Q[i][j] << "\t";
          cout << "\n\n";
        }
      }
*/

  }






  return 0;
}

double prodint(int longitud, double A[],double B[]){
  int i=0;
  double prodint=0;
  for(i=0;i<longitud;i++)
    prodint+=A[i]*B[i];
    return prodint;
}

double normalizar_vector(int longitud, double A[]){
  double vector_normalizado[longitud];
  for(int i=0;i<longitud;i++)
    vector_normalizado[i] = A[i]/sqrt(prodint(longitud,A,A));
  return 0.0;
}

void ClearScreen()
    {
    cout << string( 100, '\n' );
    }
