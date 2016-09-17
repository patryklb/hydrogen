#include <fstream>
#include <iomanip>
#include <math.h>

#include "mlib.h"
#include "dsyev.h"

using namespace std;

extern void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
double* W, double* WORK, int* LWORK, int* INFO);

static double prec = 1.0E-8;
static int max_iter = 100;
int x = 0;
int digit_prec = 4;
int n_state = 0;
double zeta = 1.0;

int n_alpha, N, res_J;
double R, R_min, R_max, R_step;
Vector alpha, d, c, e;
Matrix H, S, Q, F, V;

double* matrix_to_array(Matrix A){
    int n = A.size();
    
    double* result = new double[n * n];

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            result[i * n + j] = A[i][j];
        }
    }
    return result;
}
double* vector_to_array(Vector a){
	int n = a.size();
	double* result = new double[n];

	for(int i = 0; i < n; i++){
		result[i] = a[i];
		}
	return result;
	}
Matrix array_to_matrix(double *a, int n, int nn){
	Matrix result(n, n, 0);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < nn; j++){
			result[i][j] = a[i*n + j];
			}
		}
	return result;
	}
Vector array_to_vector(double *a, int n){
	Vector result(n);

	for(int i = 0; i < n; i++){
			result[i] = a[i];
		}
	return result;
	}
Matrix diagonalize(Matrix H, Matrix S, Vector e){
	/* Solves generalized eigenvalue problem Hv = eSv.
	 * 
	 * Input:
	 * H = hermitian matrix
	 * S = overlap matrix
	 * (H, S remain unchanged)
	 * 
	 * Output:
	 * e = vector of eigenvalues
	 * RETURN =  matrix of eigenvectors
	 * (memory allocated during the procedure) 
	 */
	 
	 /* LOCAL variables used by LAPACK routine dsyev */
	 int info;												// information about success/failure of dsyev
	 int n = H.size();										// dimension of matrices
	 int lwork = 3 * n;										// dimension of workspace for lapack routine
	 static double threshold = 1.0e-10;						// when matrix element is considered to be 0

	 
	 /* Creating arrays for dsyev */
	 double* h = matrix_to_array(H);
	 double* v = new double[n * n];
	 double* eigen_values = new double[n];
	 
	 double* work = new double [lwork];

	/* Copy S into an auxiliary matrix (dsyev destroys the matrix) */
	double* aux = new double[n*n];
	for (int i = 0; i < n; ++i ) {
		for(int j = 0; j < n; j++){
			aux[i*n + j] = S[i][j];
		}
	}
	 
	 /*  Diagonalize S  */

	dsyev_ ("V", "U", &n, aux, &n, eigen_values, work, &lwork, &info ) ;

	if (info !=0){
		fprintf (stderr, "S-matrix diagonalization failed\n");
		exit (1);
	}

  /* Keep only linearly independent combinations (within a given threshold) */
  
  int nn = 0;
  for (int  i = 0; i < n; ++i ) {
    /*  i runs on all eigenstates
       nn runs on eigenstate with nonzero eigenvalue */
    if ( eigen_values[i] > threshold) {
      for (int j = 0; j < n; ++j ) {
	aux[j+nn*n] = aux[j+i*n] / sqrt(eigen_values[i]);
      }
      ++nn;
    }
  }

  if ( nn < n ) {
     fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
   } 
  /*       Trasform H using the "aux" matrix
           V(i,j) = \sum_{k=1}^{n} H(i,k) aux(k,j),  i=1,n, j=1,nn
   */
  
  for (int i = 0; i < n; ++i ) {
    for (int j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  (int k = 0; k < n; ++k ) {
	v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
      }
    }
  }
   /* h1(i,j) = \sum_{k=1}^{n} aux(k,i) v(k,j),  i=1,nn, j=1,nn
      H' = transpose(aux)*H*aux
   */
   
  double* uno = new double[nn*nn];
  for (int i = 0; i < nn; ++i ) {
    for (int j = 0; j < nn; ++j ) {
      uno[i+j*nn] = 0.0;
      for  (int k = 0; k < n; ++k ) {
	uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
      }
    }
  }
  
  /*  Diagonalize transformed H  */
  
  dsyev_ ("V", "U", &nn, uno, &nn, eigen_values, work, &lwork, &info );

  if ( info !=0 ) {
    fprintf (stderr, "H-matrix diagonalization failed\n");
    exit (1);
  }

  /*  Back-transform eigenvectors  */

  for (int i = 0; i < n; ++i ) {
    for (int j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  (int k = 0; k < n; ++k ) {
	v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
      }
    }
  }
  

  free(uno); free(aux); free(work);
    return array_to_matrix(v, n, nn);
}
void load_data(){
    ifstream in("h2.in");
    in >> n_alpha;
    N = 2 * n_alpha;

    alpha = Vector(N);

    for(int i = 0; i < n_alpha; i++){
		in >> alpha[i];
		cout << "Gaussian #" << i << " coefficient = " << alpha[i];
		cout << endl;
	}
	in.close();

	cout << "H2 distance (a.u.): [min, max, step] >> ";
	cin >> R_min >> R_max >> R_step;
	//cout << "Energy level n >> ";
	//cin >> n_state;
	n_state = 0;
	cout << "Resonance number J >> ";
	cin >> res_J;
	  
}
void initialize_SHQ(){
    double alpha_sum, alpha_sum_ij, alpha_sum_lm, alpha_coefficient,
           nuclei_distance, R_ij, R_lm, R_ijlm;

    // initializing matrices S and H
    for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				alpha_sum = alpha[i] + alpha[j];
				alpha_coefficient = alpha[i] * alpha[j] / alpha_sum;
				nuclei_distance = pow(d[i] - d[j], 2);

                // Filling matrices S, H with terms s_ij and t_ij
                S[i][j] = pow(M_PI / alpha_sum, 1.5) * exp(-alpha_coefficient * nuclei_distance);
                H[i][j] = alpha_coefficient * S[i][j] * (6.0 - 4.0 * alpha_coefficient * nuclei_distance);

                // Now we add Coulomb repulsion term to H
                for(int nucleus = -1; nucleus < 2; nucleus += 2){
                    R_ij = abs((alpha[i] * d[i] + alpha[j] * d[j]) / alpha_sum - nucleus * R / 2.0);

                    // should we add (times 2) to both?
                    if(R_ij > prec){
                        H[i][j] = H[i][j] - 2.0 * S[i][j] * erf(sqrt(alpha_sum) * R_ij) / R_ij;

                        //if(i == 4 && j == 4 && nucleus == -1) cout << endl << R << "   " << H[i][j];
                    }
                    else{
                        H[i][j] +=  - 4.0 * M_PI / alpha_sum;
                    }


                }
		}
    }

    // initializing interaction matrix Q
    for(int i = 0; i < N; i++){
        for(int l = 0; l < N; l++){
            for(int j = 0; j < N; j++){
                for(int m = 0; m < N; m++){
                    alpha_sum_ij = alpha[i] + alpha[j];
                    alpha_sum_lm = alpha[l] + alpha[m];
                    alpha_sum = alpha_sum_ij + alpha_sum_lm;
                    alpha_coefficient = alpha_sum_ij * alpha_sum_lm / alpha_sum;

                    R_ij = abs((alpha[i] * d[i] + alpha[j] * d[j]) / alpha_sum);
                    R_lm = abs((alpha[l] * d[l] + alpha[m] * d[m]) / alpha_sum);
                    R_ijlm = R_ij - R_lm;

                    // or swap N <-> l and j<-> m
                    if(R_ijlm < prec)
						Q[N*i+l][N*j+m] = 4.0 * S[i][j] * S[l][m] / sqrt(M_PI) * sqrt(alpha_coefficient);
                    else
						Q[N*i+l][N*j+m] = 2.0 * S[i][j] * S[l][m] * erf(sqrt(alpha_coefficient) * R_ijlm) / R_ijlm;
                }
            }
        }
    }
}
void fill_fock_matrix(int mes, double R){
        // Fill Fock matrix
         for(int i = 0; i < N; i++){
			for(int l = 0; l < N; l++){
                F[i][l] = H[i][l];
                    for(int j = 0; j < N; j++){
                        for(int m = 0; m < N; m++){
                            F[i][l] += (2 * Q[N*i+j][N*l+m] - Q[N*i+l][N*j+m]) * c[j] * c[m];
                        }
                    }
			}
         }
}
int main(){
	system("clear");
	
	double new_energy, old_energy, energy_difference, nuclei_energy, zero_energy, res_energy;
	ofstream out("h2.out");

    load_data();

/*************************************************************************
Declare global workspace for matrices and vectors. Guess the initial orbital
state.
*************************************************************************/
    c = Vector(N);        // stores expansion coefficients
	d = Vector(N);        // stores position of the nuclei along z-axis
	e = Vector(N);

    F = Matrix(N, N, 0);                // Fock matrix
	S = Matrix(N, N, 0);                // overlap matrix
	V = Matrix(N, N, 0);
	H = Matrix(N, N, 0);                // full-energy matrix

	Q = Matrix(N * N, N * N, 0);  // interaction matrix

	Matrix F_temp(N, N, 0);
	Matrix S_temp(N, N, 0);

 	c[0] = 2.0;                                 // initial guess
 	c[1] = 4.0;
 	c[2] = 5.0;
 	c[3] = 5.0;
 	c[4] = 1.0;
 	c[5] = -1.0;
/***************************************************************************
Main loop over the nuclear H-H distance
***************************************************************************/
    for(R = R_min; R < R_max; R += R_step){
        //cout << endl << R << "  ";
        // setting up the nuclear positions (on the z-axis)
        for (int i = 0; i < n_alpha; i++){
			alpha[n_alpha + i] = alpha[i];
			d[i] = - R / 2.0;
			d[i + n_alpha] = R / 2.0;
		}

        initialize_SHQ();

         // Self-consistency iteration for the whole potential
         new_energy = 0.0;
		 energy_difference = 10.0;

         for(int iter = 0; (iter < max_iter && energy_difference > prec); iter++){
            fill_fock_matrix(iter, R);
          
           //diagonalize F
           // e - eigenvalues vector
           // V - matrix of eigenvectors
           // S - overlap matrixf

           F_temp = F;
           S_temp = S;
           
           //e = solve_eigen_generalized(F_temp, S_temp); // evals from highest to lowest
           //cout.precision(2);
            
           Matrix x = diagonalize(F, S, e);

          for(int i = 0; i < N; i++){
			  c[i] = x[n_state][i]; 	//maybe choose different column
			  }
          
            old_energy = new_energy;
            new_energy = 0.0;

            for(int i = 0; i < N; i++){
                for(int l = 0; l < N; l++){
                    new_energy += 2.0 * H[i][l] * c[i] * c[l];
                    for(int j = 0; j < N; j++){
                        for(int m = 0; m < N; m++){
                            new_energy += (2.0 * Q[N*i+j][N*l+m] - Q[N*i+l][N*j+m]) *
                                           c[i]*c[l]*c[j]*c[m];
                        }
                    }
                }
            }
            
            energy_difference = abs(new_energy - old_energy);
			
			cout.width(5);
			setprecision(digit_prec);
			cout.fill( ' ' );
            if(energy_difference < prec){
               // cout << "Convergence achieved [stopping]" << endl;
                if(R == R_min)
                    cout <<" [R]          [electron]      [nuclear]     [total (Ry)]   [binding energy eV]" << endl;

                if(R > prec){
                    nuclei_energy = zeta * 2.0 / R;
                    res_energy = res_J * (res_J + 1) / (2 * 1E3 * R * R);
				}
                else{
                    nuclei_energy = 0.0; 
				}
                zero_energy = -2.0;
	
                cout << setprecision(digit_prec)  << R << "\t\t" << setprecision(digit_prec)  <<  0.5 * new_energy 
					 << setprecision(digit_prec)  << "\t\t" <<  0.5 * nuclei_energy << setprecision(digit_prec) 
					 << "\t\t" <<  0.5 * (new_energy + nuclei_energy) << "\t\t" << setprecision(digit_prec)  << 0.5 * (new_energy + nuclei_energy - zero_energy) * 13.6058 << endl;
                out << R << " " << 0.5 * new_energy << " " << 0.5 * nuclei_energy << " " << (new_energy + nuclei_energy-0.15)*0.5 << " " <<
                       0.5 * (new_energy + nuclei_energy - zero_energy) * 13.6058 << endl;
            }
        }
    }
    out.close();

    return 0;
}
