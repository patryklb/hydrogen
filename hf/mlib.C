#include "mlib.h"


/**************************Class Vector********************************/
Vector::Vector(int dim) {
    v = new double [this->dim = dim];
    for (int i = 0; i < dim; i++) v[i] = 0;
}
Vector::Vector(double v0, double v1) {
    v = new double [dim = 2];
    v[0] = v0;
    v[1] = v1;
}
Vector::Vector(double v0, double v1, double v2) {
    v = new double [dim = 3];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
}
Vector::Vector(const Vector& dv) {
    v = new double [dim = dv.dim];
    for (int i = 0; i < dim; i++) v[i] = dv.v[i];
}
void Vector::resize(const int dimension) {
    double *old = new double[dim];
    for (int i = 0; i < dim; i++) old[i] = v[i];
    delete [] v;
    v = new double [dimension];
    for (int i = 0; i < dimension; i++)
        v[i] = i < dim ? old[i] : 0;
    dim = dimension;
    delete [] old;
}
void Vector::push_back(const double d) {
    resize(++dim);
    v[dim-1] = d;
}
void Vector::set(const double d0, ...) {
    va_list args;
    v[0] = d0;
    va_start(args, d0);
    for (int i = 1; i < dim; i++)
        v[i] = va_arg(args, double);
    va_end(args);
}
Vector& Vector::operator = (const Vector& dv) {
    if (this != &dv) {
        if (dim != dv.dim) {
            delete [] v;
            v = new double [dim = dv.dim];
        }
        for (int i = 0; i < dim; i++) v[i] = dv[i];
    }
    return *this;
}
Vector& Vector::operator += (const Vector& dv) {
    for (int i = 0; i < dim; i++) v[i] += dv[i];
    return *this;
}
Vector& Vector::operator -= (const Vector& dv) {
    for (int i = 0; i < dim; i++) v[i] -= dv[i];
    return *this;
}
Vector& Vector::operator *= (double d) {
    for (int i = 0; i < dim; i++) v[i] *= d;
    return *this;
}
Vector& Vector::operator /= (double d) {
    for (int i = 0; i < dim; i++) v[i] /= d;
    return *this;
}
Vector operator - (const Vector& dv) {
    int dim = dv.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = -dv[i];
    return temp;
}
Vector operator * (const Vector& dv, double d) {
    int dim = dv.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = dv[i] * d;
    return temp;
}
Vector operator * (double d, const Vector& dv) {
    int dim = dv.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = dv[i] * d;
    return temp;
}
Vector operator / (const Vector& dv, double d) {
    int dim = dv.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = dv[i] / d;
    return temp;
}
Vector operator + (const Vector& v1, const Vector& v2) {
    int dim = v1.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = v1[i] + v2[i];
    return temp;
}
Vector operator - (const Vector& v1, const Vector& v2) {
    int dim = v1.dimension();
    Vector temp(dim);
    for (int i = 0; i < dim; i++) temp[i] = v1[i] - v2[i];
    return temp;
}
double Vector::abs() {
    return std::sqrt(norm());
}

double Vector::norm() {
    double sum = 0;
    for (int i = 0; i < dim; i++) sum += v[i] * v[i];
    return sum;
}
double Vector::dot(const Vector& dv) {
    double sum = 0;
    for (int i = 0; i < dim; i++) sum += v[i] * dv[i];
    return sum;
}
std::ostream& operator<<(std::ostream& os, const Vector& dv) {
    for (int i = 0; i < dv.dim; i++) {
        os << dv.v[i];
        if (i < dv.dim-1)
            os << '\t';
        else
            os << '\n';
    }
    return os;
}
double dot(const Vector& v1, const Vector& v2) {
    double sum = 0;
    for (int i = 0; i < v1.size(); i++)
        sum += v1[i] * v2[i];
    return sum;
}

/**************************Class Matrix********************************/
static void error(const std::string str) {
    std::cerr << "cpl::Matrix error: " << str << std::endl;
    std::exit(EXIT_FAILURE);
}
static void error(const std::string str, int i) {
    std::ostringstream os;
    os << str << " " << i;
    error(os.str());
}
Matrix::Matrix(int rows, int cols, double d) {
    if (rows <= 0)
        error("bad number of rows", rows);
    if (cols <= 0)
        error("bad number of columns", cols);
    this->rows = rows;
    this->cols = cols;
    m = new double* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new double [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = d;
    }
}
Matrix::Matrix(const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;

    m = new double* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new double [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = mat[i][j];
    }
}
Matrix::~Matrix() {
    for (int i = 0; i < rows; i++)
        delete [] m[i];
    delete [] m;
}
const double* Matrix::operator[](const int row) const {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}
double* Matrix::operator[](const int row) {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}
Matrix& Matrix::operator=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols) {
            for (int i = 0; i < rows; i++)
                delete [] m[i];
            delete [] m;
            rows = mat.rows;
            cols = mat.cols;
            m = new double* [rows];
            for (int i = 0; i < rows; i++)
                m[i] = new double [cols];
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] = mat[i][j];
    }
    return *this;
}
Matrix& Matrix::operator=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] = d;
    return *this;
}
Matrix& Matrix::operator+=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols)
            error("matrix dimension mismatch");
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] += mat[i][j];
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols)
            error("matrix dimension mismatch");
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] -= mat[i][j];
    }
    return *this;
}
Matrix& Matrix::operator*=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] *= d;
    return *this;
}
Matrix& Matrix::operator/=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] /= d;
    return *this;
}
Matrix operator - (const Matrix& mat) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = -mat[i][j];
    return temp;
}
Matrix operator * (const Matrix& mat, double d) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = mat[i][j] * d;
    return temp;
}
Matrix operator * (double d, const Matrix& mat) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = d * mat[i][j];
    return temp;
}
Matrix operator / (const Matrix& mat, double d) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = mat[i][j] / d;
    return temp;
}
Matrix operator + (const Matrix& m1, const Matrix& m2) {
    int rows = m1.numRows();
    int cols = m1.numCols();
    if (rows != m2.numRows() || cols != m2.numCols())
        error("matrix dimension mismatch");
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = m1[i][j] + m2[i][j];
    return temp;
}
Matrix operator - (const Matrix& m1, const Matrix& m2) {
    int rows = m1.numRows();
    int cols = m1.numCols();
    if (rows != m2.numRows() || cols != m2.numCols())
        error("matrix dimension mismatch");
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = m1[i][j] - m2[i][j];
    return temp;
}
Matrix operator * (const Matrix& m1, const Matrix& m2) {
    int n = m1.numCols();
    if (n != m2.numRows())
        error("matrix dimension mismatch");
    int rows = m1.numRows();
    int cols = m2.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            for (int k = 0; k < n; k++)
                temp[i][j] += m1[i][k] * m2[k][j];
    return temp;
}
Matrix Matrix::transpose() {
    Matrix temp(cols, rows);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[j][i] = m[i][j];
    return temp;
}
std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            os << '\t' << mat.m[i][j];
        }
        os << '\n';
    }
    return os;
}

/**********************Numerical methods*******************************/
double pythag(double a, double b) {
	// pythagorean theorem optimized calculation
    double absa = std::abs(a);
    double absb = std::abs(b);
    if (absa > absb) {
        double ratio = absb / absa;
        return absa * std::sqrt(1 + ratio * ratio);
    } else {
        if (absb == 0.0)
            return 0.0;
        else {
            double ratio = absa / absb;
            return absb * std::sqrt(1 + ratio * ratio);
        }
    }
}
Matrix matrix_Jacobi(int d, int p, int q, double c){
	Matrix J(d, d, 0);
	double s = sqrt(1-c*c);
	for(int i = 0; i < d; i++){
		for(int j = 0; j < d; j++){
			if(i == j) J[i][j] = 1; 
			if((i == p && j == p) || (i == q && j == q)) J[i][j] = c; 
			else if(i == p && j == q) J[i][j] = s;
			else if(i == q && j == p) J[i][j] = -s;
		}
	}
	
	return J;
}
double sign(double a, double b) {
    if (b >= 0.0) {
        if (a >= 0.0)
            return a;
        else
            return -a;
    } else {
        if (a >= 0.0)
            return -a;
        else
            return a;
    }
}
Vector solve_eigen_symmetric(Matrix& A)
{   // solves eigenvalue problem A x = e x
    // returns vector of eigenvalues and replaces A by matrix of eigenvectors
    int n = A.numRows();
    Vector d(n), e(n);
    reduce_Householder(A, d, e);
    solve_TQLI(d, e, A);
    sort_eigen(d, A);
    return d;
}
void solve_eigen_generalized (Matrix& H, Matrix& S, Vector& E, Matrix& C) {
	// solves eigenvalue problem H C = E S x
    // replaces E with vector of eigenvalues and replaces H by matrix of eigenvectors
    int N = E.size();
    Matrix V = S;
    Vector s = solve_eigen_symmetric(V);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            V[i][j] /= sqrt(s[j]);
    C = V.transpose() * H * V;
    E = solve_eigen_symmetric(C);
    C = V * C;
}



void reduce_Householder(Matrix& A, Vector& d, Vector& e) {

    int n = d.dimension();

    for (int i = n - 1; i > 0; i--) {
        int l = i - 1;
        double h = 0;
        double scale = 0;
        if (l > 0) {
            for (int k = 0; k <= l; k++)
                scale += std::abs(A[i][k]);
            if (scale == 0.0)
                e[i] = A[i][l];
            else {
                for (int k = 0; k <= l; k++) {
                    A[i][k] /= scale;
                    h += A[i][k] * A[i][k];
                }
                double f = A[i][l];
                double g = (f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                A[i][l] = f - g;
                f = 0.0;
                for (int j = 0; j <= l; j++) {
                    A[j][i] = A[i][j] / h;
                    g = 0.0;
                    for (int k = 0; k <= j; k++)
                       g += A[j][k] * A[i][k];
                    for (int k = j + 1; k <= l; k++)
                       g += A[k][j] * A[i][k];
                    e[j] = g / h;
                    f += e[j] * A[i][j];
                }
                double hh = f / (h + h);
                for (int j = 0; j <= l; j++) {
                    f = A[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; k++)
                        A[j][k] -= f * e[k] + g * A[i][k];
                }
            }
        } else
            e[i] = A[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0]=0.0;
    for (int i = 0; i < n; i++) {
        if (d[i] != 0.0) {
            for (int j = 0; j < i; j++) {
                double g = 0;
                for (int k = 0; k < i; k++)
                    g += A[i][k] * A[k][j];
                for (int k = 0; k < i; k++)
                    A[k][j] -= g * A[k][i];
            }
        }
        d[i] = A[i][i];
        A[i][i] = 1.0;
        for (int j = 0; j < i; j++)
            A[j][i] = A[i][j] = 0.0;
    }
}
void solve_TQLI(Vector& d, Vector& e, Matrix& z) {

    int n = d.dimension();
    for (int i = 1; i < n; i++)
        e[i-1] = e[i];
    e[n-1] = 0.0;
    for (int l = 0; l < n; l++) {
        int iter = 0, m;
        do {
            for (m = l ; m < n-1; m++) {
                double dd = std::abs(d[m]) + std::abs(d[m+1]);
                if ((std::abs(e[m]) + dd) == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 30)
                     error("Too many iterations in solve_TQLI");
                double g = (d[l+1] - d[l]) / (2.0 * e[l]);
                double r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + sign(r, g));
                double s = 1.0;
                double c = 1.0;
                double p = 0.0;
                int i;
                for (i = m-1; i >= l; i--) {
                    double f = s * e[i];
                    double b = c * e[i];
                    e[i+1] = r = pythag(f, g);
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i+1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i+1] = g + (p = s * r);
                    g = c * r - b;
                    for (int k = 0; k < n; k++) {
                        f = z[k][i+1];
                        z[k][i+1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}
void sort_eigen(Vector& d, Matrix& z) {

    // sorts eigenvalues and eigenvector in descending order
    int n = d.dimension();
    if (z.numRows() != n || z.numCols() != n)
        error("Bad vector, matrix dimensions in sort_eigen");
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        double p = d[k];
        for (int j = i; j < n; j++)
            if (d[j] >= p)
                p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
                p = z[j][i];
                z[j][i] = z[j][k];
                z[j][k] = p;
            }
        }
    }
}


