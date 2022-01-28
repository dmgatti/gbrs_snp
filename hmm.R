################################################################################
# GBRS-SNP DO haplotype reconstruction HMM.
#
# Danielt Gatti
# dan.gatti@jax.org
# 2021-12-06
################################################################################

##### LIBRARIES #####

library(Rcpp)

##### FUNCTIONS #####

################################################################################
#  Given two arguments on a log scale, add them on a non-logged scale.
#  Arguments: x: double, value on a log scale.
#             y: double, value on a log scale.
# Return: double, the value log(exp(x) + exp(y)).
################################################################################
cppFunction("double addlog(double x, double y) {
  double retval = 0.0;
  if(!std::isfinite(x) && x < 0) {
    retval = y;
  } else if(!std::isfinite(y) && y < 0) { 
    retval =  x;
  } else if (x >= y) {
    retval = x + log1p(exp(y - x));
  } else {
    retval = y + log1p(exp(x - y));
  } /* else */

  return retval;
}")

# Notation
# alpha: forward pass values.
# beta: backward pass values.
# a: state transition probabilities.
# b: observation probabilities.
# c: scaling factor for alpha & beta to avoid underflow.
# obs: genotype observations.
# pi: initial state probabilities.

# Forward algorithm.
################################################################################
# Run the scaled Rabiner forward algorithm on a set of observations.
# pi: double[num_states] containing the initial state probabilities. Must sum to 1.0.
# a: double[num_states][num_states][num_markers] containing the state transition 
#    probabilities for each marker.
# b: double[num_genotypes][num_markers] containing the emission probabilities
#    for each genotype at each marker.
# obs: double[num_markers] containing the observed genotype at each marker for 
#      one sample.
# num_samples: int containing the number of samples.
# num_states: int containing the number of states.
# num_markers: int containing the number of markers.
################################################################################ pi: double array of length num_states containing intial state probabilities.
# a: double array with num_states rows and num_states columns containing the
#    transition probabilities.
# b: double array with num_alleles rows, num_states columns and num_markers
#    slices containing the emission probabilities P(obs|state) for each marker.
#    Each state must sum to 1.
# obs: int array of length num_markers containing the observations coded as
#      integers.
# c: double array with scaling factor.
# alpha: double array with num_genotypes row and num_markers columns containing
#        the forward values.
# num_samples: Number of samples.
# num_states: Number of genotype states.
# num_markers: Number of markers.
# num_alleles: Number of observed alleles.
# Return: void. But alpha and c will be populated with values.
cppFunction("void forward(double* pi, double** a, double** b, int* obs, double* c, double** alpha, int num_samples, int num_states, int num_markers, int num_alleles) {

  int s = 0;  /* Sample counter. */
  int g1 = 0;  /* Genotype state counter. */
  int g2 = 0;  /* Genotype state counter. */
  int m = 0;  /* Marker counter. */

  /* Initialize first marker. */
  for(g1 = 0; g1 < num_states; g1++) {
    a[g1][0] = pi[g] * b[g1][obs[0]][0];
    c[0] += a[g1][0];
  } // for(g1)
  
  /* Scale a[0]. */
  c[0] = 1.0 / c[0];
  for(g1 = 0; g1 < num_states; g1++) {
    a[g1][0] *= c[0];
  } // for(g1)
  
  for(m = 1; m < num_markers; m++) {
    c[m] = 0;
    for(g1 = 0; g1 < num_states; g1++) {
      a[g1][m] = 0;
      for(g2 = 0; g2 < num_states; g2++) {
        alpha[g1][m] += alpha[g1][m-1] * a[g2][g1]; 
      } // for(g2)
      alpha[g1][m] *= b[obs[m]][g1][m];
      c[m] += alpha[g1][m];  
    } // for(g)
    
    c[m] = 1.0 / c[m];
    for(g1 = 0; g1 < num_states; g1++) {
      alpha[g1][m] *= c[m];
    } // for(g1)
  } // for(m)
}")

cppFunction("double* backward(double* pi, double** a, double** b, int* obs, double* c, double** alpha, double** beta, int num_samples, int num_states, int num_markers, int num_alleles) {

  int s = 0;  /* Sample counter. */
  int g1 = 0;  /* Genotype state counter. */
  int g2 = 0;  /* Genotype state counter. */
  int m = 0;  /* Marker counter. */

  /* Copy scaled alpha values to beta.*/
  for(g1 = 0; g1 < num_states; g1++) {
    beta[g1][m-1] = alpha[g1][m-1];
  } // for(g1)
  
  for(m = num_markers - 2; m >= 0; m--) {
    for(g1 = 0; g1 < num_states; g1++) {
      beta[g1][m] = 0;
      for(g2 = 0; g2 < num_states; g2++) {
        beta[g1][m] += a[g1][g2] * b[g2][obs[m+1]][m+1] * beta[g2][m+1];
      } //for(g2)
      beta[g1][m] *= c[m];
    } //for(g1)
  } //for(m)
}")

##### VARIABLES #####

# Genotype states.
states = outer(LETTERS[1:8], LETTERS[1:8], paste0))
states = states[upper.tri(states, diag = TRUE]

# Founder genotypes.
founder_geno = read.csv()

# Sample genotypes.
geno = read.csv()

# Marker map.
map = read.csv()

# Number of chromosomes.
num_chr = length(map)

# Initial 36 state genotype probabilities.
init_pr = rep(1 / length(states), length(states))

# Emission probabilities.
emis_probs = vector('list', num_chr)

# Transition probabilities.
trans_pr = vector('list', num_chr)

##### PROGRAM #####

