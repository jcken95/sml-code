// stan program for hetGP
// useful for when the noise depends on the input
// e.g. sin(x) + eps_1 + eps_2*x has noise that increases with |x|

functions {
	
	matrix exp_kern( matrix Diffs, real sq_sigma, real invTheta, int N ) {
		matrix[N, N] mat;
		int ii = 0;
		mat = rep_matrix(0, N, N);
		
		for(i in 1:N){

			ii = ii + 1;

			for(j in 1:ii){

				mat[i, j] = (-1)*invTheta*pow(Diffs[i,j], 2);

				}

			}
			
			mat = mat + mat' ;
			mat = sq_sigma*exp(mat) ;
			
			return(mat);
		}

}

data {
	
	int<lower = 1> m_p;		// number regression functions for the mean
	int<lower = 1> v_p;		// number regression functions for the log-var
	int<lower = 1> N;		// number data points
	int<lower = 1> K;		// dimension of input space
	matrix[N, K] x;			// input data (should be studentised)
	matrix[N, m_p] m_H;		// design matrix (i.e. H = h(x)) (mean)
	matrix[N, v_p] v_H;		// design (log-var)
	vector[N] y;			// code outputs (noisy) - for now assume no replication ...
	matrix[N, N] Diffs;
	//vector<lower = 1>[N] a;		// vector of replication level (we can have different levels of replication here, unlike pairs of homGP)
	// leave replication out of it for now ...

	// now describe the prior


	// first the prior for the mean GP
	vector[m_p] m_beta_m;
	vector[m_p] m_beta_s; // mean and sd of mean function parameters for the mean


	real<lower = 0> m_a_theta;

	real<lower = 0> m_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> m_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> m_b_sigma;

	real<lower = 0> m_nugget; // useful for stabilising the inversion of matrices

	// next the prior for the log-variance GP
	vector[v_p] v_beta_m;
	vector[v_p] v_beta_s; // mean and sd of mean function parameters for the mean

	real<lower = 0> v_a_theta;
	real<lower = 0> v_b_theta; //correlation lengthscales (will use the form 1/theta^2)
	
	real<lower = 0> v_a_sigma; //params of dist for sigma^2 (log-var)
	real<lower = 0> v_b_sigma;

	real<lower = 0> v_nugget_a; // quantifies the noise of the noise
	real<lower = 0> v_nugget_b;
	// we might set the nugget terms to be 10^(-4)
	
}


parameters {

vector[m_p] m_beta;		// parameters of mean function for mean
vector[v_p] v_beta;		// parameters of mean function for variance

real<lower = 0> m_sigma;	// mean scale param
real<lower = 0> v_sigma;	// log var scale param

real<lower = 0> m_theta; 	// length scale parameters
real<lower = 0> v_theta;

real<lower = 0> v_nugget;

vector[N] logLambda;

}

model {

vector[N] mu;			// the "actual" mean of the top-level GP
vector[N] v_mu;	// mean of the latent log-variance

matrix[N, N] m_var;		// variance matrix of the mean
matrix[N, N] v_var;		// variance matrix of the log-variance

matrix[N, N] m_var_chol;	// cholesky them bois
matrix[N, N] v_var_chol;

vector[N] lambda;		// latent variance


// first produce the variance

for(i in 1:N){
	mu[i] = m_H[i, ] * m_beta;
	v_mu[i] = v_H[i, ] * v_beta;
}

v_var = exp_kern(Diffs, square(v_sigma), pow(v_theta, -2), N) + diag_matrix(rep_vector(v_nugget*v_nugget, N));



lambda = exp(logLambda);

m_var = exp_kern(Diffs, square(m_sigma), pow(m_theta, -2), N) + diag_matrix(lambda) + diag_matrix(rep_vector(m_nugget*m_nugget, N));

y ~ multi_normal(mu, m_var); // top level statement

// prior beliefs
print(logLambda);

logLambda ~ multi_normal(v_mu, v_var); // statement about the log-variance

for(i in 1:m_p){
	m_beta[i] ~ normal(m_beta_m[i], m_beta_s[i]);
}

for(i in 1:v_p){
	v_beta[i] ~ normal(v_beta_m[i], v_beta_s[i]);
}

m_theta ~ gamma(m_a_theta, m_b_theta);
v_theta ~ gamma(v_a_theta, v_b_theta);

m_sigma ~ inv_gamma(m_a_sigma, m_b_sigma);
v_sigma ~ inv_gamma(v_a_sigma, v_b_sigma);

v_nugget ~ inv_gamma(v_nugget_a, v_nugget_b);

}


