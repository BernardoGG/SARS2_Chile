// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n_obs;
  real<lower=0> weeks[n_obs];
  int<lower=0> all_lins[n_obs];
  int<lower=0> pers_lins[n_obs];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0, upper=1> theta;
  real lambda;
}

transformed parameters {
  real<lower=0, upper=1> prop_lins[n_obs];
  for(i in 1:n_obs){
      prop_lins[i] = theta*exp(-lambda*weeks[i]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
    pers_lins ~ binomial(all_lins, prop_lins);
}

