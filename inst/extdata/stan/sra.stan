functions {

  // This function is to facilitate solving the system to obtain lambda_max
  vector solve_lambda(vector y, vector theta, real[] x_r, int[] x_i) {
    real A = theta[1];
    real S = theta[2];
    //real lambda_max = exp(y[1]);
    real lambda_max = y[1];
    vector[1] z;
    z[1] = exp(pow(A + (S / (lambda_max - S)), -1)) - lambda_max;
    return z;
  }
  
  // return the Euclidean norm of a vector
  real vector_norm(vector x) {
	    
	    real i = 0.0;
	    
	    for (j in 1:num_elements(x))
	        i += pow(x[j], 2.0);
	        
	    return pow(i, 0.5);
  }
  
  /*
	* Return the number of non-zero elements in an
	* integer vector
	*/
	int num_nonzero(int[] y) {
		int np = 0;
		for (n in 1:size(y))
			if (y[n] > 0)
				np += 1;
		return np;
	}

} // end of functions


data {

  // Dimensions
  int<lower=1> n_species; // the number of species (71)
  int          n_month; // (12)
  int<lower=1> n_method; // the number of fishing methods (4)
  int<lower=n_method> n_fishery_group; // the number of fishery groups (19)
  int<lower=1,upper=n_species> n_species_group; // the number of species groups (30)
  int<lower=1,upper=n_species> n_cm_group; // the number of cryptic mortality groups (5)
  int<lower=1> n_fishery_group_m[n_method]; // the number of fishery groups for each method
  int<lower=1> idx_method_fg[n_fishery_group]; // index vector indcating which method applies to which group
  int<lower=1,upper=n_species_group> species_group_s[n_species]; // species groups for vulnerability

  // Cryptic mortality priors
  int<lower=1,upper=n_cm_group> cm_group_s[n_species]; // cryptic mortality group for each species
  int<lower=1> n_cm; // number of cryptic mortality parameters (48)
  vector<lower=0>[n_cm] cm_par1; // beta priors
  vector<lower=0>[n_cm] cm_par2;
  int<lower=1,upper=n_cm> cm_id[n_fishery_group,n_cm_group];

  // Biological priors
  // N_BP
  int<lower=0,upper=1> n_breeding_pairs_type[n_species];
  vector<lower=0>[n_species] n_breeding_pairs_p1;
  vector<lower=0>[n_species] n_breeding_pairs_p2;
  // P_B
  vector<lower=0>[n_species] p_breeding_p1;
  vector<lower=0>[n_species] p_breeding_p2;
  // A_curr
  vector<lower=0>[n_species] age_breeding_current_p1;
  vector<lower=0>[n_species] age_breeding_current_p2;
  // S_curr
  int<lower=0,upper=1> adult_survival_current_type[n_species];
  vector<lower=0>[n_species] adult_survival_current_p1;
  vector<lower=0>[n_species] adult_survival_current_p2;
  // S_opt
  int<lower=0,upper=1> adult_survival_opt_type[n_species];
  vector<lower=0>[n_species] adult_survival_opt_p1;
  vector<lower=0>[n_species] adult_survival_opt_p2;
  // p_EEZ
  matrix[n_species, n_month] p_eez;
  // p_nest
  matrix[n_species, n_month] p_nest;
  
  // Observed fishing events
  int<lower=1> n_i; // number of observed fishing events
  int<lower=1,upper=n_month> month_i[n_i];
  int<lower=1,upper=n_method> method_i[n_i]; // BLL, SLL, SN, trawl
  int<lower=1,upper=n_fishery_group> fishery_group_i[n_i];
  int<lower=1,upper=n_species> species_i[n_i];
  int<lower=1,upper=n_species_group> species_group_i[n_i];
  int<lower=1,upper=n_cm_group> cm_group_i[n_i];
  vector<lower=0>[n_i] overlap_i;
  int<lower=0> live_captures_i[n_i];
  int<lower=0> dead_captures_i[n_i];
  int<lower=0> captures_i[n_i];

  // All fishing events
  int<lower=1> n_j; // number of observed and unobserved fishing events
  int<lower=1,upper=n_month> month_j[n_j];
  int<lower=1,upper=n_method> method_j[n_j]; // BLL, SLL, SN, trawl
  int<lower=1,upper=n_fishery_group> fishery_group_j[n_j];
  int<lower=1,upper=n_species> species_j[n_j];
  int<lower=1,upper=n_species_group> species_group_j[n_j];
  int<lower=1,upper=n_cm_group> cm_group_j[n_j];
  vector<lower=0>[n_j] overlap_j;
  
  real<lower=1> n_years;

  real<lower=0,upper=1> psi;

} // end of data


transformed data {

  real x_r[0];
  int x_i[0];
  vector[1] y_guess;
  y_guess[1] = 1.0;

} // end of transformed data


parameters {

  vector[n_method] log_q_intercept;
  vector[n_fishery_group-n_method] log_q_method_raw;
  vector[n_species_group-1] log_q_species_raw;
  
  real<lower=0> tau;
  matrix[n_fishery_group,n_species] eps_gs;

  vector[n_fishery_group] beta_live_capture_fg;
  vector[n_species_group] beta_live_capture_sg;
  
  vector<lower=0,upper=1>[n_species-sum(n_breeding_pairs_type)] n_breeding_pairs_raw1_s; // log-uniform
  vector<lower=0>[sum(n_breeding_pairs_type)] n_breeding_pairs_raw2_s; // log-normal
  vector[n_species] p_breeding_raw_s; // logit-normal
  
} // end of parameters


transformed parameters {

  // catchability (probability of
  // observing a capture)
  vector<lower=0>[n_fishery_group] q_sg[n_species];
  
  // catchability predictands
  vector[n_species_group] log_q_z;
  vector[n_fishery_group] log_q_g;
  
  // probabilities of live capture
  vector<lower=0,upper=1>[n_fishery_group] p_live_capture_zg[n_species_group];
  vector<lower=0,upper=1>[n_i] p_live_captures_i;

  // biological parameters
  vector<lower=0>[n_species] n_breeding_pairs_s;
  vector<lower=0>[n_species] n_adults_s;
  vector<lower=0>[n_month]   n_vuln_adults_s[n_species];
  vector<lower=0,upper=1>[n_species] p_breeding_s;
  
  // expected captures
  vector[n_i] mu_observed_captures_i;
  vector[n_i] mu_live_captures_i;
  vector[n_i] mu_dead_captures_i;

  // Species group vulnerability expressed in log-space 
  // and forced to sum to zero
  {
    vector[n_fishery_group] q_method;
 
    // log-catchability per species group
    // with sum-to-zero constraint
    for (z in 1:(n_species_group-1)) {
      log_q_z[z] = log_q_species_raw[z];
    }
    log_q_z[n_species_group] = -sum(log_q_species_raw);
    
    // log-catchability per fishery group
    // with sum-to-zero constraint
    // within each method
    {
      int c1 = 1;
      int c2 = 1;
      for(mm in 1:n_method){

        int c3 = c2 + n_fishery_group_m[mm] - 1;

        if(n_fishery_group_m[mm] > 1) {
          real tmp_sum = 0;
          for(ii in 0:(n_fishery_group_m[mm] - 2)) {
            q_method[c2+ii] = log_q_method_raw[c1 + ii];
            tmp_sum += log_q_method_raw[c1 + ii];
          } 
          q_method[c3] = -tmp_sum;
        } else {
          q_method[c3] = 0;
        }
        c1 += n_fishery_group_m[mm] - 1;
        c2 += n_fishery_group_m[mm];
      }
    }
    log_q_g = log_q_intercept[idx_method_fg] + q_method;
    
    // back transform catchability onto natural scale
    for (g in 1:n_fishery_group) {
      for (s in 1:n_species) {
        q_sg[s,g] = inv_logit(log_q_z[species_group_s[s]] + log_q_g[g] + eps_gs[g,s]);
      }
    }
    
    // probability of live capture as a function of coefficients
    for (g in 1:n_fishery_group) {
      for (z in 1:n_species_group) {
        p_live_capture_zg[z,g] = inv_logit(beta_live_capture_fg[g] + beta_live_capture_sg[z]);
      }
    }
  }
  
  // predicted biological parameters
  {
    int c1 = 1;
    int c2 = 1;

    for (s in 1:n_species) {
      real L;
      real U;
      // N_BP: 0=log-uniform or 1=lognormal
      if (n_breeding_pairs_type[s] == 0) {
        L = log(n_breeding_pairs_p1[s]);
	    U = log(n_breeding_pairs_p2[s]);
        n_breeding_pairs_s[s] = exp(L + (U - L) * n_breeding_pairs_raw1_s[c1]);
        c1 = c1 + 1;
      } else {
        n_breeding_pairs_s[s] = n_breeding_pairs_raw2_s[c2];
        c2 = c2 + 1;
      }
      // P_B: logit-normal
      p_breeding_s[s] = inv_logit(p_breeding_raw_s[s]);
      
      // N_adults (per annum)
      n_adults_s[s] = 2.0 * n_breeding_pairs_s[s] / p_breeding_s[s];
      
      // N_vuln_adults (per species and month)
      for(m in 1:n_month){
        n_vuln_adults_s[s, m] = n_adults_s[s] * p_eez[s, m] * (1 - p_breeding_s[s] * p_nest[s, m]);
      }
    }
  }

  {
      real density_overlap_i;
      
      for (i in 1:n_i) {
          
        density_overlap_i = fmax(1e-8, overlap_i[i] * n_vuln_adults_s[species_i[i], month_i[i]]);
        
        mu_observed_captures_i[i] = q_sg[species_i[i],fishery_group_i[i]] * density_overlap_i;
        
        p_live_captures_i[i]  = p_live_capture_zg[species_group_i[i], fishery_group_i[i]];
        
        mu_live_captures_i[i] = mu_observed_captures_i[i] * p_live_captures_i[i];
        mu_dead_captures_i[i] = mu_observed_captures_i[i] * (1.0 - p_live_captures_i[i]);
      }
  }

} // end of transformed parameters


model {

  // priors on catchability parameters
  log_q_intercept   ~ normal(0.0, 10.0);
  log_q_method_raw  ~ normal(0.0, 2.0);
  log_q_species_raw ~ normal(0.0, 2.0);

  // random effects error
  // on catchability
  for (g in 1:n_fishery_group) {
    eps_gs[g] ~ normal(0.0, tau);
  }
  tau ~ cauchy(0.0, 1.0);
  
  // priors on estimated biological 
  // parameters
  {
      int c1 = 1;
      real mu;
      real sigma;
      
      for (s in 1:n_species) {

        // N_BP: 0=log-uniform or 1=lognormal
        if (n_breeding_pairs_type[s] == 1) {
          mu = n_breeding_pairs_p1[s];
          sigma = n_breeding_pairs_p2[s];
          n_breeding_pairs_raw2_s[c1] ~ lognormal(log(mu), sigma);
          c1 = c1 + 1;
        }
        
        // P_B - logit-normal
        mu = p_breeding_p1[s];
        sigma = p_breeding_p2[s];
        p_breeding_raw_s[s] ~ normal(logit(mu), sigma / (mu * (1.0 - mu)));
      }
  }

  // priors on live capture coefficients
  beta_live_capture_fg ~ logistic(0.0, 1.0);
  beta_live_capture_sg ~ logistic(0.0, 1.0);

  // likelihood for total observed captures
  captures_i ~ poisson(mu_observed_captures_i);
  
  // likelihood for conditional live capture
  // given capture
  {
      int n_nz = num_nonzero(captures_i);
      
      int captures_k[n_nz];
      int live_captures_k[n_nz];
      
      real p_live_captures_k[n_nz];
      
      int k = 1;
      
      for (i in 1:n_i) {
          if (captures_i[i] > 0) {
              live_captures_k[k]   = live_captures_i[i];
              captures_k[k]        = captures_i[i];
              p_live_captures_k[k] = p_live_captures_i[i];
              k = k + 1;
          }
      }
      
      live_captures_k ~ binomial(captures_k, p_live_captures_k);
  }

} // end of model


generated quantities {

  // summary traces
  vector[7] traces;
  
  // priors
  real prior_tau;
  vector[n_fishery_group] prior_beta_live_capture_fg;
  vector[n_species_group] prior_beta_live_capture_sg;
  matrix[n_species_group,n_fishery_group] prior_p_live_capture_zg;
  vector[n_species] prior_n_breeding_pairs_s;
  vector[n_species] prior_p_breeding_s;
  vector[n_species] prior_n_adults_s;
  matrix[n_species,n_fishery_group] prior_p_live_capture_sg;

  // posterior predictive checking
  real p_survive_capture;
  vector[n_i] pred_live_captures_i;
  vector[n_i] pred_dead_captures_i;
  
  // cryptic mortality
  vector<lower=0,upper=1>[n_cm] p_observable;
  
  // back calculation of vulnerability
  // (probability of capture)
  matrix[n_species,n_fishery_group] vulnerability_sg;

  // Captures and deaths
  vector[n_species]       observed_captures_s;
  vector[n_method]        observed_captures_m;
  vector[n_fishery_group] observed_captures_g;
  vector[n_method]        observed_captures_sm[n_species];
  vector[n_fishery_group] observed_captures_sg[n_species];
  vector[n_species]       observed_deaths_s;
  vector[n_method]        observed_deaths_m;
  vector[n_fishery_group] observed_deaths_g;
  vector[n_method]        observed_deaths_sm[n_species];
  vector[n_fishery_group] observed_deaths_sg[n_species];
  
  vector[n_species]       captures_s;
  vector[n_method]        captures_m;
  vector[n_fishery_group] captures_g;
  vector[n_method]        captures_sm[n_species];
  vector[n_fishery_group] captures_sg[n_species];
  vector[n_species]       deaths_s;
  vector[n_method]        deaths_m;
  vector[n_fishery_group] deaths_g;
  vector[n_method]        deaths_sm[n_species];
  vector[n_fishery_group] deaths_sg[n_species];
  
  vector[n_species]       cryptic_deaths_s;
  vector[n_method]        cryptic_deaths_m;
  vector[n_fishery_group] cryptic_deaths_g;
  vector[n_method]        cryptic_deaths_sm[n_species];
  vector[n_fishery_group] cryptic_deaths_sg[n_species];
  
  // Mortality check
  vector[n_species] mortality_in_bounds_s;

  // Outputs for Risk Atlas - these outputs need to be standardised across all models being used by Risk Atlas
  matrix[n_species,n_fishery_group] p_observable_sg;
  matrix[n_species,n_fishery_group] p_live_capture_sg;
  matrix[n_species,n_fishery_group] p_survive_capture_sg;
  
  // rmax, PST, risk ratio
  vector[n_species] age_breeding_current_s;
  vector[n_species] adult_survival_current_s;
  vector[n_species] adult_survival_opt_s;
  vector[n_species] r_max_s;
  vector[n_species] pst_s;
  
  vector[n_species]                 risk_ratio_s;
  vector[n_method]                  risk_ratio_m;
  vector[n_fishery_group]           risk_ratio_g;
  matrix[n_species,n_method]        risk_ratio_sm;
  matrix[n_species,n_fishery_group] risk_ratio_sg;
  
  vector[n_species]                 q_risk_s;
  vector[n_method]                  q_risk_m;
  vector[n_fishery_group]           q_risk_g;
  
  // Cryptic multiplier
  for (i in 1:n_cm) {
    p_observable[i] = beta_rng(cm_par1[i], cm_par2[i]);
  }
  
  // calculate vulnerability
  for(s in 1:n_species){
    for(g in 1:n_fishery_group){
      vulnerability_sg[s,g] = q_sg[s,g] / p_observable[cm_id[g, cm_group_s[s]]];
    }
  }
  
  // Prior checking and simulated parameters
  prior_tau = cauchy_rng(0.0, 1.0);
  for (g in 1:n_fishery_group) {
    prior_beta_live_capture_fg[g] = logistic_rng(0.0, 1.0);
  }
  for (s in 1:n_species_group) {
    prior_beta_live_capture_sg[s] = logistic_rng(0.0, 1.0);
  }
  for (g in 1:n_fishery_group) {
    for (z in 1:n_species_group) {
      prior_p_live_capture_zg[z,g] = inv_logit(prior_beta_live_capture_fg[g] + prior_beta_live_capture_sg[z]);
    }
  }
  
  {
      real mu;
      real sigma;
      vector[2] theta;
      vector[1] lambda_max;
        
      for (s in 1:n_species) {
        
        // N_BP - log-uniform or log-normal
        mu = n_breeding_pairs_p1[s];
        sigma = n_breeding_pairs_p2[s];
        if (n_breeding_pairs_type[s] == 0) {
          prior_n_breeding_pairs_s[s] = exp(uniform_rng(log(mu), log(sigma)));
        } else {
          prior_n_breeding_pairs_s[s] = lognormal_rng(log(mu), sigma);
        }
        
        // P_B - logit normal
        mu = p_breeding_p1[s];
        sigma = p_breeding_p2[s];
        prior_p_breeding_s[s] = inv_logit(normal_rng(logit(mu), sigma / (mu * (1.0 - mu))));
        
        // N_adults
        prior_n_adults_s[s] = 2.0 * n_breeding_pairs_s[s] / prior_p_breeding_s[s];
        
        // S_curr - uniform or logit-normal
        mu = adult_survival_current_p1[s];
        sigma = adult_survival_current_p2[s];
        if (adult_survival_current_type[s] == 0) {
          adult_survival_current_s[s] = uniform_rng(mu, sigma);
        } else {
          adult_survival_current_s[s] = inv_logit(normal_rng(logit(mu), sigma / (mu * (1.0 - mu))));
        }
        
        // S_opt - uniform or logit-normal
        mu = adult_survival_opt_p1[s];
        sigma = adult_survival_opt_p2[s];
        if (adult_survival_opt_type[s] == 0) {
          adult_survival_opt_s[s] = uniform_rng(mu, sigma);
        } else {
          adult_survival_opt_s[s] = inv_logit(normal_rng(logit(mu), sigma / (mu * (1.0 - mu))));
        }
        
        // A_curr - uniform
        age_breeding_current_s[s] = uniform_rng(age_breeding_current_p1[s], age_breeding_current_p2[s]);

        // not sure what this is
        p_survive_capture = uniform_rng(0, 1);

        // PST
        theta[1] = age_breeding_current_s[s];
        theta[2] = adult_survival_opt_s[s];
        lambda_max = algebra_solver(solve_lambda, y_guess, theta, x_r, x_i);
        r_max_s[s] = log(lambda_max[1]);
        pst_s[s] = 0.5 * psi * r_max_s[s] * n_adults_s[s];
      }
  }
  
  // Outputs as matrices by species and fishery group for Risk Atlas
  for (s in 1:n_species) {
    for (g in 1:n_fishery_group) {
    
      prior_p_live_capture_sg[s,g] = prior_p_live_capture_zg[species_group_s[s],g];

      p_observable_sg[s,g] = p_observable[cm_id[g,cm_group_s[s]]];
      p_live_capture_sg[s,g] = p_live_capture_zg[species_group_s[s],g];
      p_survive_capture_sg[s,g] = p_survive_capture;
    }
  }
    
  // initialise captures and deaths
  for (s in 1:n_species) {      
    for (g in 1:n_fishery_group) {
        
      observed_captures_sg[s,g] = 0.0;
      observed_deaths_sg[s,g]   = 0.0;
        
      captures_sg[s,g] = 0.0;
      deaths_sg[s,g]   = 0.0;
      
      cryptic_deaths_sg[s,g] = 0.0;
    }
    for (m in 1:n_method) {
        
      observed_captures_sm[s,m] = 0.0;
      observed_deaths_sm[s,m]   = 0.0;
        
      captures_sm[s,m] = 0.0;
      deaths_sm[s,m]   = 0.0;
      
      cryptic_deaths_sm[s,m] = 0.0;
    }
  }

  // posterior predictive checking
  for (i in 1:n_i) {
    
    pred_live_captures_i[i] = poisson_rng(mu_live_captures_i[i]);
    pred_dead_captures_i[i] = poisson_rng(mu_dead_captures_i[i]);
  }

  // get annual captures and deaths
  {
    real density_overlap_j;
    real p_live_captures_j;
    real mu_observed_captures_j;
    real mu_observed_deaths_j;
    real mu_captures_j;
    real mu_deaths_j;
    real captures_j;
    real deaths_j;
    
      for (j in 1:n_j) {
        
        density_overlap_j = overlap_j[j] * n_vuln_adults_s[species_j[j], month_j[j]];
        
        p_live_captures_j = p_live_capture_zg[species_group_j[j], fishery_group_j[j]];
        
        // expected observed captures per year
        mu_observed_captures_j = q_sg[species_j[j], fishery_group_j[j]] * density_overlap_j / n_years;
        
        // expected observed deaths per year
        mu_observed_deaths_j = mu_observed_captures_j * (1.0 - p_live_captures_j);
        
        // simulate obsered captures and deaths
        if (mu_observed_captures_j > 0) {
            
            // posterior prediction of average total captures and deaths
            captures_j = poisson_rng(n_years * mu_observed_captures_j) / n_years;
            deaths_j   = poisson_rng(n_years * mu_observed_deaths_j) / n_years;
            
            // sum captures
            observed_captures_sm[species_j[j], method_j[j]] += captures_j;
            observed_captures_sg[species_j[j], fishery_group_j[j]] += captures_j;
            
            // sum deaths
            observed_deaths_sm[species_j[j], method_j[j]] += deaths_j;
            observed_deaths_sg[species_j[j], fishery_group_j[j]] += deaths_j;
        }
        
        // expected total captures per year
        mu_captures_j = vulnerability_sg[species_j[j], fishery_group_j[j]] * density_overlap_j / n_years;
        
        // expected total deaths per year (observed + unobserved)
        mu_deaths_j = mu_captures_j * (1.0 - p_live_captures_j * p_survive_capture);
        
        // simulate total captures and deaths
        if (mu_captures_j > 0) {
            
            // posterior prediction of average total captures and deaths
            captures_j = poisson_rng(n_years * mu_captures_j) / n_years;
            deaths_j   = poisson_rng(n_years * mu_deaths_j) / n_years;
            
            // sum captures
            captures_sm[species_j[j], method_j[j]] += captures_j;
            captures_sg[species_j[j], fishery_group_j[j]] += captures_j;
            
            // sum deaths
            deaths_sm[species_j[j], method_j[j]] += deaths_j;
            deaths_sg[species_j[j], fishery_group_j[j]] += deaths_j;
            
            // expected cryptic deaths 
            cryptic_deaths_sm[species_j[j], method_j[j]]        += deaths_j * 
                                                                   (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            cryptic_deaths_sg[species_j[j], fishery_group_j[j]] += deaths_j * 
                                                                   (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
        }
      }
  }
  
  // summmary results by species
  for (s in 1:n_species) {
      
    observed_captures_s[s] = sum(observed_captures_sm[s]);
    observed_deaths_s[s]   = sum(observed_deaths_sm[s]);
    
    captures_s[s]       = sum(captures_sm[s]);
    deaths_s[s]         = sum(deaths_sm[s]);
    
    cryptic_deaths_s[s] = sum(cryptic_deaths_sm[s]);
    
    risk_ratio_s[s]     = deaths_s[s] / pst_s[s];
    
    // riks of capture per unit effort
    q_risk_s[s]         = inv_logit(log_q_z[species_group_s[s]]);
    
    for (m in 1:n_method) {
      risk_ratio_sm[s,m] = deaths_sm[s,m] / pst_s[s];
    }
    for (g in 1:n_fishery_group) {
      risk_ratio_sg[s,g] = deaths_sg[s,g] / pst_s[s];
    }
  }
  
  // summary results by method
  for (m in 1:n_method) {
      
    observed_captures_m[m] = sum(observed_captures_sm[,m]);
    observed_deaths_m[m]   = sum(observed_deaths_sm[,m]);
    
    captures_m[m]       = sum(captures_sm[,m]);
    deaths_m[m]         = sum(deaths_sm[,m]);
    
    cryptic_deaths_m[m] = sum(cryptic_deaths_sm[,m]);
    
    // riks of capture per unit effort
    q_risk_m[m]         = inv_logit(log_q_intercept[m]);
    
    // mean risk accross species weighted by species-specific deaths
    risk_ratio_m[m] = sum((to_vector(risk_ratio_sm[,m]) .* deaths_s) / sum(deaths_s));
  }
  
  // summary results by fishery_group
  for (g in 1:n_fishery_group) {
      
    observed_captures_g[g] = sum(observed_captures_sg[,g]);
    observed_deaths_g[g]   = sum(observed_deaths_sg[,g]);
    
    captures_g[g]       = sum(captures_sg[,g]);
    deaths_g[g]         = sum(deaths_sg[,g]);
    
    cryptic_deaths_g[g] = sum(cryptic_deaths_sg[,g]);
   
    // riks of capture per unit effort
    q_risk_g[g]         = inv_logit(log_q_g[g]);
    
    // mean risk accross species weighted by species-specific deaths
    risk_ratio_g[g] = sum((to_vector(risk_ratio_sg[,g]) .* deaths_s) / sum(deaths_s));
  }
  
  // Mortality constraint
  {
      real total_deaths_s;
      
      for (s in 1:n_species) {
          
        total_deaths_s = (1.0 - adult_survival_current_s[s]) * n_adults_s[s];
        
        mortality_in_bounds_s[s] = deaths_s[s] >= total_deaths_s ? 0 : 1;
      }
  }
  
  // summary trace outputs
  traces[1] = vector_norm(log_q_intercept);
  traces[2] = vector_norm(log_q_method_raw);
  traces[3] = vector_norm(log_q_species_raw);
  traces[4] = tau;
  traces[5] = vector_norm(append_row(beta_live_capture_fg, beta_live_capture_sg));
  traces[6] = vector_norm(append_row(n_breeding_pairs_raw1_s, n_breeding_pairs_raw2_s));
  traces[7] = vector_norm(p_breeding_raw_s);
  
} // end of generated quantities
