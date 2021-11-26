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
  int<lower=n_method>                n_fishery_group; // the number of fishery groups (19)
  int<lower=1,upper=n_species>       n_species_group; // the number of species groups (30)
  int<lower=1>                       n_fishery_group_m[n_method]; // the number of fishery groups for each method
  int<lower=1>                       idx_method_fg[n_fishery_group]; // index vector indcating which method applies to which group
  
  int<lower=1,upper=n_species_group> species_group_s[n_species]; // species groups for vulnerability

  // Cryptic mortality multipliers (log-normal)
  vector[n_species] cm_longline_par1;
  vector[n_species] cm_longline_par2;
  vector[n_species] cm_warp_par1;
  vector[n_species] cm_warp_par2;
  vector[n_species] cm_net_par1;
  vector[n_species] cm_net_par2;
  
  int n_capture_type; // 3
  int n_capture_type_group; // 5
  
  int<lower=1,upper=n_capture_type_group> capture_type_group_s[n_species];

  // Biological priors
  // N_BP
  int<lower=0,upper=2> n_breeding_pairs_type[n_species];
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
  vector<lower=0>[n_i] overlap_i;
  int<lower=0> live_captures_i[n_i];
  int<lower=0> dead_captures_i[n_i];
  int<lower=0> net_captures_i[n_i];
  int<lower=0> warp_captures_i[n_i];
  int<lower=0> other_captures_i[n_i];
  int<lower=0> captures_i[n_i];
  
  // Species capture type groups
  int capture_type_group_i[n_i];

  // All fishing events
  int<lower=1> n_j; // number of observed and unobserved fishing events
  int<lower=1,upper=n_month> month_j[n_j];
  int<lower=1,upper=n_method> method_j[n_j]; // BLL, SLL, SN, trawl
  int<lower=1,upper=n_fishery_group> fishery_group_j[n_j];
  int<lower=1,upper=n_species> species_j[n_j];
  int<lower=1,upper=n_species_group> species_group_j[n_j];
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
  
  vector<lower=0>[n_species] n_breeding_pairs_raw_s; 
  vector[n_species] p_breeding_raw_s; // logit-normal
  
  simplex[n_capture_type] p_capture_type_t[n_capture_type_group];
  
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
        q_sg[s,g] = exp(log_q_z[species_group_s[s]] + log_q_g[g] + eps_gs[g,s]);
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
    real L;
    real U;
      
    for (s in 1:n_species) {
      
      // N_BP: 0=log-uniform or 1=lognormal
      if (n_breeding_pairs_type[s] == 0) {
        L = log(n_breeding_pairs_p1[s]);
	    U = log(n_breeding_pairs_p2[s]);
        n_breeding_pairs_s[s] = exp(L + (U - L) * n_breeding_pairs_raw_s[s]);
      } else {
        n_breeding_pairs_s[s] = n_breeding_pairs_raw_s[s];
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
          n_breeding_pairs_raw_s[c1] ~ lognormal(log(mu), sigma);
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
  for (i in 1:n_i) {
      if (captures_i[i] > 0) {
          live_captures_i[i] ~ binomial(captures_i[i], p_live_captures_i[i]);
      }
  }
  
  // likelihood for trawl capture types
  for (i in 1:n_i) {
      if (method_i[i] == 4) {
          
          int captures_method_i = net_captures_i[i] + warp_captures_i[i] + other_captures_i[i];
          
          if (captures_method_i > 0) {
          
            net_captures_i[i]   ~ binomial(captures_method_i, p_capture_type_t[capture_type_group_i[i], 1]);
            warp_captures_i[i]  ~ binomial(captures_method_i, p_capture_type_t[capture_type_group_i[i], 2]);
            other_captures_i[i] ~ binomial(captures_method_i, p_capture_type_t[capture_type_group_i[i], 3]);
          }
      }
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
  
  vector<lower=0,upper=1>[n_fishery_group] p_live_capture_sg[n_species];

  // posterior predictive checking
  real p_survive_capture;
  vector[n_i] observed_captures_i;
  vector[n_i] observed_live_captures_i;
  vector[n_i] observed_dead_captures_i;
  
  // back calculation of vulnerability
  // (probability of capture)
  matrix[n_species, n_fishery_group] vulnerability_sg;
  // (probability of observation)
  matrix<lower=0,upper=1>[n_species, n_fishery_group] p_observable_sg;

  // Captures and deaths
  vector[n_species]       observed_captures_s = rep_vector(0.0, n_species);
  vector[n_method]        observed_captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observed_captures_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       observed_deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        observed_deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observed_deaths_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_species]       observable_captures_s = rep_vector(0.0, n_species);
  vector[n_method]        observable_captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observable_captures_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       observable_deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        observable_deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observable_deaths_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_species]       captures_s = rep_vector(0.0, n_species);
  vector[n_method]        captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] captures_g = rep_vector(0.0, n_fishery_group);
  vector[n_method]        deaths_sm[n_species] = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_fishery_group] deaths_sg[n_species] = rep_array(rep_vector(0.0, n_fishery_group), n_species);
  vector[n_species]       deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] deaths_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_species]       cryptic_deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        cryptic_deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] cryptic_deaths_g = rep_vector(0.0, n_fishery_group);
  
  // Mortality check
  vector[n_species] mortality_in_bounds_s;
  
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
  for(s in 1:n_species){
    for(g in 1:n_fishery_group){
        prior_p_live_capture_sg[s,g] = prior_p_live_capture_zg[species_group_s[s],g];
    }
  }
  
  // calculate vulnerability
  {
      real cm;
      
      for(s in 1:n_species){
        for(g in 1:n_fishery_group){
            
            p_live_capture_sg[s,g] = p_live_capture_zg[species_group_s[s],g];
            
            // Longline
            if (idx_method_fg[g] == 1 || idx_method_fg[g] == 2) {
                
                cm = lognormal_rng(log(cm_longline_par1[s]), cm_longline_par2[s]);
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg[s,g] + (1 - p_live_capture_sg[s,g]) * cm);
                
            } else { 
            
            // Set net
            if (idx_method_fg[g] == 3) {
                
                cm = 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg[s,g] + (1 - p_live_capture_sg[s,g]) * cm);
                
            } else { 
            
            // Trawl
            if (idx_method_fg[g] == 4) {
                
                vulnerability_sg[s,g] = q_sg[s,g] * p_live_capture_sg[s,g];
                
                // Net
                cm = 1.0 + lognormal_rng(log(cm_net_par1[s]), cm_net_par2[s]);
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg[s,g]) * p_capture_type_t[capture_type_group_s[s], 1] * cm;
                
                // Warp
                cm = lognormal_rng(log(cm_warp_par1[s]), cm_warp_par2[s]);
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg[s,g]) * p_capture_type_t[capture_type_group_s[s], 2] * cm;
                
                // Other
                cm = 1.0;
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg[s,g]) * p_capture_type_t[capture_type_group_s[s], 3] * cm;
            }}}
            
            p_observable_sg[s,g] = q_sg[s,g] / vulnerability_sg[s,g];
        }
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

  // get total captures and deaths
  // from observed overlap
  for (i in 1:n_i) {
    
    observed_captures_i[i]      = poisson_rng(mu_observed_captures_i[i]);
    observed_live_captures_i[i] = poisson_rng(mu_live_captures_i[i]);
    observed_dead_captures_i[i] = poisson_rng(mu_dead_captures_i[i]);
    
    // sum captures
    observed_captures_s[species_i[i]]       += observed_captures_i[i];
    observed_captures_g[fishery_group_i[i]] += observed_captures_i[i];
    observed_captures_m[method_i[i]]        += observed_captures_i[i];
    
    // sum captures
    observed_deaths_s[species_i[i]]       += observed_dead_captures_i[i];
    observed_deaths_g[fishery_group_i[i]] += observed_dead_captures_i[i];
    observed_deaths_m[method_i[i]]        += observed_dead_captures_i[i];
  }

  // get annual captures and deaths
  // from total overlap
  {
    real density_overlap_j;
    real p_live_captures_j;
    real mu_observable_captures_j;
    real mu_observable_deaths_j;
    real mu_captures_j;
    real mu_deaths_j;
    real captures_j;
    real deaths_j;
    
      for (j in 1:n_j) {
        
        density_overlap_j = overlap_j[j] * n_vuln_adults_s[species_j[j], month_j[j]];
        
        p_live_captures_j = p_live_capture_zg[species_group_j[j], fishery_group_j[j]];
        
        // expected observable captures per year
        mu_observable_captures_j = q_sg[species_j[j], fishery_group_j[j]] * density_overlap_j / n_years;
        
        // expected observable deaths per year
        mu_observable_deaths_j = mu_observable_captures_j * (1.0 - p_live_captures_j);
        
        // simulate obsered captures and deaths
        if (mu_observable_captures_j > 0) {
            
            // posterior prediction of average total captures and deaths
            captures_j = poisson_rng(n_years * mu_observable_captures_j) / n_years;
            deaths_j   = poisson_rng(n_years * mu_observable_deaths_j) / n_years;
            
            // sum captures
            observable_captures_s[species_j[j]]       += captures_j;
            observable_captures_g[fishery_group_j[j]] += captures_j;
            observable_captures_m[method_j[j]]        += captures_j;
            
            // sum deaths
            observable_deaths_s[species_j[j]]       += deaths_j;
            observable_deaths_g[fishery_group_j[j]] += deaths_j;
            observable_deaths_m[method_j[j]]        += deaths_j;
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
            captures_s[species_j[j]]       += captures_j;
            captures_g[fishery_group_j[j]] += captures_j;
            captures_m[method_j[j]]        += captures_j;
            
            // sum deaths
            deaths_sm[species_j[j], method_j[j]]        += deaths_j;
            deaths_sg[species_j[j], fishery_group_j[j]] += deaths_j;
            
            deaths_s[species_j[j]]       += deaths_j;
            deaths_g[fishery_group_j[j]] += deaths_j;
            deaths_m[method_j[j]]        += deaths_j;
            
            // expected cryptic deaths 
            cryptic_deaths_s[species_j[j]]       += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            cryptic_deaths_g[fishery_group_j[j]] += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            cryptic_deaths_m[method_j[j]]        += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
        }
      }
  }
  
  // summmary results by species
  for (s in 1:n_species) {
    
    // risk per species
    risk_ratio_s[s]     = deaths_s[s] / pst_s[s];
    
    // riks of capture per unit effort
    q_risk_s[s] = inv_logit(mean(log_q_intercept) + log_q_z[species_group_s[s]]);
    
    for (m in 1:n_method) {
      risk_ratio_sm[s,m] = deaths_sm[s,m] / pst_s[s];
    }
    for (g in 1:n_fishery_group) {
      risk_ratio_sg[s,g] = deaths_sg[s,g] / pst_s[s];
    }
  }
  
  // summary results by method
  for (m in 1:n_method) {
    
    // riks of capture per unit effort
    q_risk_m[m]         = inv_logit(log_q_intercept[m]);
    
    // mean risk accross species weighted by species-specific deaths
    risk_ratio_m[m] = sum((to_vector(risk_ratio_sm[,m]) .* deaths_s) / sum(deaths_s));
  }
  
  // summary results by fishery_group
  for (g in 1:n_fishery_group) {
   
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
  traces[6] = vector_norm(n_breeding_pairs_raw_s);
  traces[7] = vector_norm(p_breeding_raw_s);
  
} // end of generated quantities
