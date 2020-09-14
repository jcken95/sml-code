// HetGP but with an SML level added into it
// the cheap code will simply be modelled by a homoscedastic GP 
// because I'm not massively interested in the noise here
functions {
	
	vector[,] diffs( matrix X1, matrix X2, int n1, int n2, int K ) {
		
		// X1 might be the cheap design
		// X2 might be the expensive design
		
		// or they might be the same

		// Ni is the size of design i


		vector[K] D[n1,n2];
		vector[K] tmp1;
		vector[K] tmp2;
		
		
		
		for(i in 1:n1){
			for(j in 1:n2){
				tmp1 = (X1[i,]' - X2[j,]');
				D[i,j] = tmp1 .* tmp1;
			}
		}
		
		return(D);
	}

	matrix exp_kern( vector[,] D, real sq_sigma, vector omega, int n1, int n2,    int K  ) {

		matrix[n1,n2] mat;
		vector[K] tmp1;
		real tmp2;
		
		
		int ii = 0;

		mat = rep_matrix(0, n1, n2);
		if(n1 != n2){
			for(i in 1:n1){

				//ii = ii + 1; // ii = 1:(N-1)

				for(j in 1:n2){ // just compute the lower triangle

					tmp1 = D[i,j];	

					mat[i,j] =  (-1)*dot_product( omega, tmp1 );

				}

	
			}
		}
		if(n1 == n2){
			for(i in 1:n1){

				ii = ii + 1; // ii = 1:(N-1)

				for(j in ii:n2){ // just compute the lower triangle

					tmp1 = D[i,j];	

					mat[i,j] =  (-1)*dot_product( omega, tmp1 );

				}

	
			}
	
			mat = mat + mat';

		}


		mat = sq_sigma*exp(mat);

		
		return(mat);
	}

}

data {
	
	int<lower = 1> m_p;		// number regression functions for the mean
	int<lower = 1> v_p;		// number regression functions for the log-var (expensive)
	int<lower = 1> N;		// number data points
	int<lower = 1> n_c;		// number cheap points
	int<lower = 2> n_e;		// number expensive points
	int<lower = 2> K;		// dimension of input space
	matrix[n_e, K] x_e;		// expensive input variables (should be studentised)
	matrix[n_c, K] x_c;		// expensive input variables (should be studentised)
	matrix[n_c, 2*m_p] m_H;		// design matrix (i.e. H = h(x)) of the mean, should be in block form
	matrix[n_e, v_p] v_H;		// design (log-var of expensive code)
	vector[n_c] y_c;			// code outputs (noisy) - for now assume no replication ...
	vector[n_e] y_e;
	
	// now describe the priors of GP hyperparams


	// first the prior for the mean GP
	vector[m_p] m_beta_m;
	vector[m_p] m_beta_s; // mean and sd of mean function parameters for the mean


	vector<lower = 0>[K] m_a_theta;
	vector<lower = 0>[K] m_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> m_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> m_b_sigma;

	real<lower = 0> m_nugget; // useful for stabilising the inversion of matrices

	// next the prior for the log-variance GP
	vector[v_p] v_beta_m;
	vector[v_p] v_beta_s; // mean and sd of mean function parameters for the mean

	vector<lower = 0>[K] v_a_theta;
	vector<lower = 0>[K] v_b_theta; //correlation lengthscales (will use the form 1/theta^2)
	
	real<lower = 0> v_a_sigma; //params of dist for sigma^2 (log-var)
	real<lower = 0> v_b_sigma;

	real<lower = 0> v_nugget_a; // quantifies the noise of the noise
	real<lower = 0> v_nugget_b;
	// we might set the nugget terms to be 10^(-4)

	// prior for the cheap mean GP
	
	vector[m_p] c_beta_m;
	vector[m_p] c_beta_s; // mean and sd of mean function parameters for the mean


	vector<lower = 0>[K] c_a_theta;
	vector<lower = 0>[K] c_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> c_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> c_b_sigma;

	real<lower = 0> c_nugget_a;	// quantifies the noise of the cheap mean
	real<lower = 0> c_nugget_b;	// just assume constant noise on this for simplicity	
	
	real m_rho;
	real<lower = 0> s_rho;
}
transformed data{
		
	vector[N] y; // all the code outputs concatenated
	vector[K] D_ee[n_e,n_e]; // their pairwise differences in block form
	vector[K] D_cc[n_c,n_c];
	vector[K] D_ce[n_c,n_e];

	D_ee = diffs(x_e, x_e,  n_e, n_e, K);
	D_cc = diffs(x_c, x_c, n_c, n_c, K);
	D_ce = diffs(x_c, x_e, n_c, n_e, K);
	y = append_row(y_c, y_e);	

}


parameters {

vector[m_p] m_beta;		// parameters of mean function for mean
vector[v_p] v_beta;		// parameters of mean function for variance
vector[m_p] c_beta;		// parameters of cheap mean fn

real<lower = 0> m_sigma;	// exp mean scale param
real<lower = 0> v_sigma;	// log var scale param
real<lower = 0> c_sigma;	// cheap mean scale param

vector<lower = 0>[K] m_theta; 	// length scale parameters
vector<lower = 0>[K] v_theta;
vector<lower = 0>[K] c_theta;

real<lower = 0> v_nugget;	// variance of the log variance
real<lower = 0> c_nugget;	// variance of the cheap outputs

real rho;

vector[n_e] logLambda; // estimated log-variance at the expensive design points

}

model{

vector[N] mu;
vector[n_c] c_mu;
vector[n_e] m_mu;
vector[n_e] v_mu;
matrix[N, N] Sigma_mat;
matrix[n_e, n_e] Var_Y_e;
matrix[n_c, n_c] Var_Y_c;
matrix[n_e, n_c] Cov_c_e;
matrix[n_e, n_e] v_var; // cov structure for the latent log variance process
matrix[N, N] var_cholesky;
matrix[n_e, n_e] var_var_cholesky;

// construct the means

for(i in 1:n_c){
// m_H is in block form!

	c_mu[i] = m_H[i,1:m_p]*c_beta;	
		
}
for(i in 1:n_e){
	m_mu[i] = m_H[i + n_c - n_e ,(1+m_p):2*m_p] * (m_beta + rho*c_beta); 	// m_H is in block form!
	v_mu[i] = v_H[i,] * v_beta;
}

mu = append_row(c_mu, m_mu);
// construct the variance matrix

Var_Y_c = exp_kern(D_cc, square(c_sigma), exp(-2*log(c_theta)), n_c, n_c, K);

// I can be clever and order the x_c such that the intersection with x_e is at end
Var_Y_e = rho * rho * exp_kern(D_ee, square(c_sigma), exp(-2*log(c_theta)), n_e, n_e, K) + exp_kern(D_ee, square(m_sigma), exp(-2*log(m_theta)), n_e, n_e, K) + diag_matrix(exp(logLambda)); 


Cov_c_e = rho * exp_kern(D_ce, square(c_sigma), exp(-2*log(c_theta)), n_c, n_e, K)' ;

Var_Y_c = Var_Y_c + diag_matrix(rep_vector(c_nugget, n_c));

v_var = exp_kern(D_ee, square(v_sigma), exp(-2*log(v_theta)), n_e, n_e, K) + diag_matrix(rep_vector(v_nugget, n_e));

Sigma_mat = append_row( append_col(Var_Y_c, Cov_c_e') , append_col(Cov_c_e, Var_Y_e) ); 

var_cholesky = cholesky_decompose(Sigma_mat);
print(dims(Sigma_mat));
y ~ multi_normal_cholesky(mu, var_cholesky); // statement about the observed code outputs
// the above is <kind of> two GPs

// prior beliefs

var_var_cholesky = cholesky_decompose(v_var);
logLambda ~ multi_normal_cholesky(v_mu, var_var_cholesky); // statement about the log-variance (latent structure)

for(i in 1:m_p){
	m_beta[i] ~ normal(m_beta_m[i], m_beta_s[i]);
}

for(i in 1:m_p){
	c_beta[i] ~ normal(c_beta_m[i], c_beta_s[i]);
}

for(i in 1:v_p){
	v_beta[i] ~ normal(v_beta_m[i], v_beta_s[i]);
}

for(k in 1:K){
	m_theta[k] ~ gamma(m_a_theta[k], m_b_theta[k]);
	c_theta[k] ~ gamma(c_a_theta[k], c_b_theta[k]);
	v_theta[k] ~ gamma(v_a_theta[k], v_b_theta[k]);
}

m_sigma ~ inv_gamma(m_a_sigma, m_b_sigma);
c_sigma ~ inv_gamma(c_a_sigma, c_b_sigma);
v_sigma ~ inv_gamma(v_a_sigma, v_b_sigma);

c_nugget ~ inv_gamma(c_nugget_a, c_nugget_b);
v_nugget ~ inv_gamma(v_nugget_a, v_nugget_b);

rho ~ normal(m_rho, s_rho); // a priori we expect the codes to be giving a similar output
print(rho);

}

