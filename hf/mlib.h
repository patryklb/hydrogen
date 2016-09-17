#ifndef ML_HPP_INCLUDED
#define ML_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdarg>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include <complex>
#include <cstdarg>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <limits>

/**************************Class Vector********************************/
class Vector {
  public:
    Vector(int dim = 1);
    Vector(double, double);
    Vector(double, double, double);
    Vector(const Vector&);
    ~Vector() { delete [] v; }
   
    int dimension() const { return dim; }
    int size() const { return dim; }
    void resize(const int);
    void push_back(const double);
    void set(const double, ...);
    
    const double operator[](const int i) const { return v[i]; }
    double& operator[](const int i) { return v[i]; }
    Vector& operator = (const Vector&);
    Vector& operator += (const Vector&);
    Vector& operator -= (const Vector&);
    Vector& operator *= (double);
    Vector& operator /= (double);

    double abs();
    double norm();
    double dot(const Vector& dv);

    friend std::ostream& operator<<(std::ostream& os, const Vector& dv);

  private:
    int dim;
    double *v;
};

inline Vector operator + (const Vector& dv) {
    return dv;
}
Vector operator - (const Vector&);
Vector operator * (const Vector&, double);
Vector operator * (double, const Vector&);
Vector operator / (const Vector& dv, double d);
Vector operator + (const Vector&, const Vector&);
Vector operator - (const Vector&, const Vector&);
double dot(const Vector&, const Vector&);

/***************************Class Matrix*******************************/
class Matrix {
  public:
    Matrix(int rows=1, int cols=1, double d=0);
    Matrix(const Matrix&);
    ~Matrix();
    int numRows() const { return rows; }
    int size() const { return rows; }
    int numCols() const { return cols; }
    const double* operator[](const int) const;
    double* operator[](const int);
    Matrix& operator=(const Matrix&);
    Matrix& operator=(const double);
    Matrix& operator += (const Matrix&);
    Matrix& operator -= (const Matrix&);
    Matrix& operator *= (double);
    Matrix& operator /= (double);
    Matrix transpose();

    friend std::ostream& operator<<(std::ostream&, const Matrix&);

  private:
    int rows;
    int cols;
    double **m;
};

inline Matrix operator + (const Matrix& m) {
    return m;
}

Matrix operator - (const Matrix&);
Matrix operator * (const Matrix&, double);
Matrix operator * (double, const Matrix&);
Matrix operator / (const Matrix& m, double d);
Matrix operator + (const Matrix& m1, const Matrix& m2);
Matrix operator - (const Matrix& m1, const Matrix& m2);
Matrix operator * (const Matrix& m1, const Matrix& m2);

inline double sign(double, double);
double pythag(double, double);


Matrix matrix_Jacobi(int, int, int, double);



//void solve_Gauss_Jordan(Matrix& A, Matrix& B);
//void solve_LU_decompose(Matrix& A, Matrix& B);
void reduce_Householder(Matrix&, Vector&, Vector&);
void solve_TQLI(Vector&, Vector&, Matrix&);
void sort_eigen(Vector& d, Matrix& z);
Vector solve_eigen_symmetric(Matrix& A);
void solve_eigen_generalized (Matrix& H, Matrix& S, Vector& E, Matrix& C);
//void singular_value_decompose(Matrix& a, Vector& w, Matrix& v);
//void minimize_BFGS(Vector& p, const double gtol, int& iter, double& fret,
//                  double (*func)(Vector&), void (*dfunc)(Vector&, Vector&));





#endif /* ML_HPP_INCLUDED */
