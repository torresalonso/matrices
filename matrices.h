#include <math.h>
/*Librería para trabajo con matrices.*/
typedef struct{int reng; int cols; double matriz[100][100];} matriz_mn;
typedef struct{int dimension; double matriz[100][100];} matriz_cuadrada;
typedef struct{int dimension; double componentes[100];} vector;

double prodint(vector u, vector v);
vector normalizar_vector(vector v);

matriz_mn matriz_transpuesta(matriz_mn A);
matriz_cuadrada matriz_transpuesta(matriz_cuadrada A);

matriz_mn reiniciar_matriz(matriz_mn A);
matriz_cuadrada  reiniciar_matriz(matriz_cuadrada A);
matriz_mn producto_matrices(matriz_mn A, matriz_mn B);
matriz_cuadrada producto_matrices(matriz_cuadrada A, matriz_cuadrada B);

matriz_cuadrada QR_getQ(matriz_cuadrada A, bool unitarios = 1);
matriz_cuadrada QR_getR(matriz_cuadrada A, matriz_cuadrada Q);

//Eigenvalores
matriz_cuadrada QR_eigenvalores(matriz_cuadrada A);
//---------------------------------------------Definición de funciones

double prodint(vector u, vector v){
  int i=0;
  double prodint=0;
  for(i=0;i<u.dimension;i++)
    prodint+=u.componentes[i]*v.componentes[i];
    return prodint;
}

vector normalizar_vector(vector v){
  vector vector_normalizado;
  vector_normalizado.dimension = v.dimension;
  for(int i = 0; i < v.dimension; i++)
    vector_normalizado.componentes[i] = v.componentes[i]/sqrt(prodint(v,v));
  return vector_normalizado;
}

matriz_mn matriz_transpuesta(matriz_mn A){
  int i = 0, j = 0;
  matriz_mn transpuesta;
  transpuesta.reng = A.cols;
  transpuesta.cols = A.reng;

  for(i = 0; i < transpuesta.reng; i++)
    for(j = 0; j < transpuesta.cols; j++)
      transpuesta.matriz[i][j] = A.matriz[j][i];

  return transpuesta;
}

matriz_cuadrada matriz_transpuesta(matriz_cuadrada A){
  int i = 0, j = 0;
  matriz_cuadrada transpuesta;
  transpuesta.dimension = A.dimension;

  for(i = 0; i < transpuesta.dimension; i++)
    for(j = 0; j < transpuesta.dimension; j++)
      transpuesta.matriz[i][j] = A.matriz[j][i];

  return transpuesta;
}

//Producto de matrices
matriz_mn reiniciar_matriz(matriz_mn A){
  for (int i = 0; i < A.reng; i++)
    for (int j = 0; j < A.cols; j++)
      A.matriz[i][j] = 0;
  return A;

}
matriz_cuadrada  reiniciar_matriz(matriz_cuadrada A){
  for (int i = 0; i < A.dimension; i++)
    for (int j = 0; j < A.dimension; j++)
      A.matriz[i][j] = 0;
  return A;
}

matriz_mn producto_matrices(matriz_mn A, matriz_mn B){
  int i = 0, j = 0, k = 0;
  matriz_mn producto;
  producto.reng = A.reng;
  producto.cols = A.cols;

  producto = reiniciar_matriz(producto);

  for(i = 0; i < producto.reng; i++)
    for(j = 0; j < producto.cols; j++){
      //copiar el renglón i-ésimo de A
      for(k = 0; k < producto.cols; k++)
        producto.matriz[i][j] += A.matriz[i][k]*B.matriz[k][j];
    }
  return producto;
}

matriz_cuadrada producto_matrices(matriz_cuadrada A, matriz_cuadrada B){
  int i = 0, j = 0, k = 0;
  matriz_cuadrada producto;
  producto.dimension = A.dimension;
  producto = reiniciar_matriz(producto);

  for(i = 0; i < producto.dimension; i++)
    for(j = 0; j < producto.dimension; j++){
      //copiar el renglón i-ésimo de A
      for(k = 0; k < producto.dimension; k++)
        producto.matriz[i][j] += A.matriz[i][k]*B.matriz[k][j];
    }

  return producto;
}

//Matrices ortogonales.
matriz_cuadrada QR_getQ(matriz_cuadrada A, bool unitarios){
  int n=0, i=0,j=0, k=0, l=0, flag=0;
  //---------------gram-schmidt---------------
  //iniciar la matriz ortogonal y el vector auxiliar
  matriz_cuadrada Q;
  vector aux, a, q;

  n = Q.dimension = aux.dimension = a.dimension = q.dimension = A.dimension;

  for(i = 0; i < n; i++){
    aux.componentes[i] = 0;
    for (j = 0; j < n; j++)
      Q.matriz[i][j] = 0;
  }

  //q1 = a1
  for(i = 0; i < n; i++)
    Q.matriz[i][0] = q.componentes[i] = A.matriz[i][0];
  if(unitarios)
    for(i = 0; i < n; i++)
      Q.matriz[i][0] = A.matriz[i][0]/sqrt(prodint(q,q));

  //--Obtener los vectores q_j, para j>1:
  /*el vector q almacena las columnas ya ortogonalizadas*/
  for(j = 1; j < n; j++){

    //tomar la j-ésima columna de la matriz A
    for(i = 0; i < n; i++)
      a.componentes[i] = A.matriz[i][j];

    //reinicia el vector auxiliar
    for(i = 0; i < n; i++)
      aux.componentes[i] = 0;

    //tomar las columnas desde la 0 hasta la (j-1)-ésima de la matriz ortogonal
    for(k = 0; k < j; k++){

      //copiar el vector q
      for(i = 0; i < n; i++)
        q.componentes[i]=Q.matriz[i][k];

      //el vector aux guardará la combinacion lineal
      for(i = 0; i < n; i++)
        aux.componentes[i] += ((prodint(a,q))/prodint(q,q))*q.componentes[i];
    }//al cabo de este for, ya tenemos la combinación lineal de los vectores q

    //Ortogonalizamos la j-ésima columna de A:
    for(i = 0; i < n; i++)
      Q.matriz[i][j] = q.componentes[i] = a.componentes[i]-aux.componentes[i];

    //Normalizamos el vector q
    if(unitarios)
      for(i = 0; i < n; i++)
        Q.matriz[i][j] /= sqrt(prodint(q,q));


  }
  return Q;
}

matriz_cuadrada QR_getR(matriz_cuadrada A, matriz_cuadrada Q){
  matriz_cuadrada R;
  R.dimension = A.dimension;

  R = producto_matrices(matriz_transpuesta(Q), A);
  return R;
}

matriz_cuadrada QR_eigenvalores(matriz_cuadrada A, int iteraciones){
  int i=0;
  matriz_cuadrada Ak = A, Q;
  for(i=0;i<iteraciones;i++){
    Q = QR_getQ(Ak, 1);
    Ak = producto_matrices(QR_getR(Ak,Q),Q);
  }
  return Ak;
}
