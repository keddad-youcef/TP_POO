#include <vector>

//TD3//

// LocalLinearAlgebra class data & methods.

class LocalLinearAlgebra {

public:

    struct Matrix {
      Matrix (int size_i, int size_j, double val) :nb_row(size_i), nb_col(size_j), V(nb_row*nb_col, val){}
      int nb_row,nb_col;
      double init_val;
      std::vector<double> V;
      void add_value(int i, int j, double val);
    };
    using Vector =  std::vector<double>;

    
    void static mult(Matrix const &A, Vector const &b, Vector &c);
    void static copy(Vector const &a, Vector &b);
    void static axpy(int a, Vector const &x, Vector &b);
    double static norm2(Vector const &r);
};