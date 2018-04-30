/*Librería para trabajo con matrices.*/
typedef struct{int dimension; double matriz[100][100];} matriz_cuadrada;
typedef struct{int dimension; double componentes[100];} vector;

double prodint(vector u, vector v);
vector normalizar_vector(vector v);
matriz_cuadrada gram_schmidt(matriz_cuadrada A, bool unitarios = 1);

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

matriz_cuadrada gram_schmidt(matriz_cuadrada A, bool unitarios){
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
