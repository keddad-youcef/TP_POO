#include <iostream>
#include <vector>
#include <math.h>
#include "LocalLinearAlgebra.h"
#include "AlienMock.h"

//TD3//


//mult() method definition
void LocalLinearAlgebra::mult(Matrix const &A, Vector const &b, Vector &c) {
  for (int i = 0; i < A.nb_row; i++) {
    for (int j = 0; j < A.nb_col; j++)
    {
      c[j] += A.V[i * A.nb_col + j] * b[j];
    }
    
  }
}

//copy() method definition
void LocalLinearAlgebra::copy(Vector const &a, Vector &b) {
  b = a;
}

//axpy() method definition
void LocalLinearAlgebra::axpy(int a, Vector const &x, Vector &b) {
  for (int i = 0; i < x.size(); i++)
  {
    b[i] += a * x[i];
  }
  
}

//norm2() method definition
double LocalLinearAlgebra::norm2(Vector const &r) {
  double norm = 0.0;
  for (int i = 0; i < r.size(); i++)
  {
    norm += r[i] * r[i];
  }
  return sqrt(norm);
}


//add_value() method definition
void LocalLinearAlgebra::Matrix::add_value(int i,int j,double val) {
  int index = (i)*nb_col + j ;
  V[index]= val; 

}
