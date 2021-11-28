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
  int n_species;                      // number of species (71)
  int n_month;                        // number of months (12)
  int n_method;                       // number of fishing methods (4)
  int n_fishery_group;                // number of fishery groups
  int n_species_group;                // number of species groups
  int n_fishery_group_m[n_method];    // number of fishery groups for each method
  int idx_method_fg[n_fishery_group]; // index vector indcating which method applies to which group
  int n_capture_type;                 // number of trawl capture types (3)
  int n_capture_type_group;           // number of trawl capture type groups
  
  int n_breeding_pairs_type[n_species];
  int n_breeding_pairs_type_0;
  int n_breeding_pairs_type_1;
  int n_breeding_pairs_type_2;
  
  // Cryptic mortality multipliers (log-normal)
  vector[n_species] cm_longline_par1;
  vector[n_species] cm_longline_par2;
  vector[n_species] cm_warp_par1;
  vector[n_species] cm_warp_par2;
  vector[n_species] cm_net_par1;
  vector[n_species] cm_net_par2;
  
  // structural groups per species
  int<lower=1,upper=n_species_group>      species_group_s[n_species];      // species groups for vulnerability
  int<lower=1,upper=n_capture_type_group> capture_type_group_s[n_species]; // species groups for capture type

  // Biological priors
  // N_BP
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
  
  // Biological covariates
  // p_EEZ
  matrix[n_species, n_month] p_eez;
  // p_nest
  matrix[n_species, n_month] p_nest;
  
  // Observed fishing events
  int<lower=1>                            n_i; 
  int<lower=1,upper=n_month>              month_i[n_i];
  int<lower=1,upper=n_method>             method_i[n_i]; 
  int<lower=1,upper=n_fishery_group>      fishery_group_i[n_i];
  int<lower=1,upper=n_species>            species_i[n_i];
  int<lower=1,upper=n_species_group>      species_group_i[n_i];
  int<lower=1,upper=n_capture_type_group> capture_type_group_i[n_i];
  
  // Observed captures per event
  int<lower=0> live_captures_i[n_i];
  int<lower=0> dead_captures_i[n_i];
  int<lower=0> net_captures_i[n_i];
  int<lower=0> warp_captures_i[n_i];
  int<lower=0> other_captures_i[n_i];
  int<lower=0> captures_i[n_i];
  
  // Observed density overlap
  vector<lower=0>[n_i] overlap_i;
  
  // All fishing events
  int<lower=1>                       n_j; 
  int<lower=1,upper=n_month>         month_j[n_j];
  int<lower=1,upper=n_method>        method_j[n_j];
  int<lower=1,upper=n_fishery_group> fishery_group_j[n_j];
  int<lower=1,upper=n_species>       species_j[n_j];
  int<lower=1,upper=n_species_group> species_group_j[n_j];
  
  // Total density overlap
  vector<lower=0>[n_j] overlap_j;

  // Data years
  real<lower=1> n_years;

  // PST tuning parameter
  real<lower=0,upper=1> psi;

} // end of data
transformed data {

  real x_r[0];
  int x_i[0];
  vector[1] y_guess;
  y_guess[1] = 1.0;

} // end of transformed data
parameters {

  // catchability predictors
  vector[n_method]                   log_q_intercept;
  vector[n_fishery_group - n_method] log_q_method_raw;
  vector[n_species_group - 1]        log_q_species_raw;
  
  // catchability
  // random effect
  real<lower=0> tau;
  matrix[n_fishery_group,n_species_group] eps_gz;

  // Predictors of live capture
  // probability
  vector[n_fishery_group] beta_live_capture_g;
  vector[n_species_group] beta_live_capture_z;
  
  // Predictors of biological values
  vector<lower=0,upper=1>[n_breeding_pairs_type_0] n_breeding_pairs_raw_0; 
  vector<lower=0>[n_breeding_pairs_type_1] n_breeding_pairs_raw_1;
  vector<lower=0>[n_breeding_pairs_type_2] n_breeding_pairs_raw_2;
  vector[n_species] p_breeding_raw_s;
  
  // Distribution of trawl captures across capture types
  // per capture type group
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
  
  // expected observed captures
  vector<lower=0>[n_i] mu_observed_captures_i;
  vector<lower=0>[n_i] mu_live_captures_i;
  vector<lower=0>[n_i] mu_dead_captures_i;
  vector<lower=0>[n_i] mu_net_captures_i;
  vector<lower=0>[n_i] mu_warp_captures_i;

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
        q_sg[s,g] = exp(log_q_z[species_group_s[s]] + log_q_g[g] + eps_gz[g,species_group_s[s]]);
      }
    }
    
    // probability of live capture as a function of coefficients
    for (g in 1:n_fishery_group) {
      for (z in 1:n_species_group) {
        p_live_capture_zg[z,g] = inv_logit(beta_live_capture_g[g] + beta_live_capture_z[z]);
      }
    }
  }
  
  // predicted biological parameters
  {
      real parA;
      real parB;
      
      int k0 = 1;
      int k1 = 1;
      int k2 = 1;
      
      for (s in 1:n_species) {
          
        parA = n_breeding_pairs_p1[s];
        parB = n_breeding_pairs_p2[s];
      
        // N_BP:
        if (n_breeding_pairs_type[s] == 0) {
          n_breeding_pairs_s[s] = parA + (parB - parA) * n_breeding_pairs_raw_0[k0];
          k0 += 1;
        }
        if (n_breeding_pairs_type[s] == 1) {
          n_breeding_pairs_s[s] = n_breeding_pairs_raw_1[k1];
          k1 += 1;
        }
        if (n_breeding_pairs_type[s] == 2) {
          n_breeding_pairs_s[s] = n_breeding_pairs_raw_2[k2];
          k2 += 1;
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
        
        // Trawl only
        mu_net_captures_i[i]  = method_i[i] == 4 ? mu_observed_captures_i[i] * p_capture_type_t[capture_type_group_i[i], 1] : 0.0;
        mu_warp_captures_i[i] = method_i[i] == 4 ? mu_observed_captures_i[i] * p_capture_type_t[capture_type_group_i[i], 2] : 0.0;
      }
  }
} // end of transformed parameters


model {

  // priors on catchability parameters
  log_q_intercept   ~ std_normal();
  log_q_method_raw  ~ std_normal();
  log_q_species_raw ~ std_normal();

  // random effects error
  // on catchability
  //for (g in 1:n_fishery_group) {
  //  eps_gz[g] ~ normal(0.0, tau); // to_vector
  //}
  to_vector(eps_gz) ~ std_normal();
  tau ~ cauchy(0.0, 1.0);
  
  // priors on estimated biological 
  // parameters
  {
      real parA;
      real parB;
      
      int k0 = 1;
      int k1 = 1;
      int k2 = 1;
      
      for (s in 1:n_species) {
          
        parA = n_breeding_pairs_p1[s];
        parB = n_breeding_pairs_p2[s];

        // N_BP:
        if (n_breeding_pairs_type[s] == 0) {
          n_breeding_pairs_raw_0[k0] ~ uniform(0, 1);
          k0 += 1;
        }
        if (n_breeding_pairs_type[s] == 1) {
          n_breeding_pairs_raw_1[k1] ~ lognormal(log(parA) - square(parB) / 2, parB);
          k1 += 1;
        }
        if (n_breeding_pairs_type[s] == 2) {
          n_breeding_pairs_raw_2[k2] ~ normal(parA, parB);
          k2 += 1;
        }
        
        // P_B:
        parA = p_breeding_p1[s];
        parB = p_breeding_p2[s];
        p_breeding_raw_s[s] ~ normal(logit(parA), parB / (parA * (1.0 - parA)));
        // Jacobian
        //target += - log(p_breeding_s[s]) - log(1 - p_breeding_s[s]);
      }
  }

  // priors on live capture coefficients
  beta_live_capture_g ~ std_normal();
  beta_live_capture_z ~ std_normal();

  // likelihood for total observed captures
  captures_i ~ poisson(mu_observed_captures_i);
  
  // likelihood for live capture
  // given capture   
  for (i in 1:n_i) {
      if (captures_i[i] > 0) {
          live_captures_i[i] ~ binomial(captures_i[i], p_live_captures_i[i]);
      }
  }
  
  // likelihood for trawl capture types
  // given capture
  for (i in 1:n_i) {
      if (method_i[i] == 4) {
          
          if (captures_i[i] > 0) {
          
            net_captures_i[i]   ~ binomial(captures_i[i], p_capture_type_t[capture_type_group_i[i], 1]);
            warp_captures_i[i]  ~ binomial(captures_i[i], p_capture_type_t[capture_type_group_i[i], 2]);
            other_captures_i[i] ~ binomial(captures_i[i], p_capture_type_t[capture_type_group_i[i], 3]);
          }
      }
  }
} // end of model
generated quantities {

  // summary traces
  vector[7] traces;
  
  // Priors for biological parameters
  vector[n_species] prior_n_breeding_pairs_raw_s; 
  vector[n_species] prior_p_breeding_raw_s;
  
  // priors for derived biological parameters
  vector<lower=0>[n_species] prior_n_breeding_pairs_s;
  vector<lower=0>[n_species] prior_n_adults_s;
  vector<lower=0>[n_month]   prior_n_vuln_adults_s[n_species];
  vector<lower=0,upper=1>[n_species] prior_p_breeding_s;

  // posterior predictive checking
  vector[n_i] observed_captures_i;
  vector[n_i] observed_live_captures_i;
  vector[n_i] observed_dead_captures_i;
  vector[n_i] observed_net_captures_i;
  vector[n_i] observed_warp_captures_i;
  
  // back calculation of vulnerability
  // (probability of capture)
  matrix[n_species, n_fishery_group] vulnerability_sg;
  // (probability of observation)
  matrix[n_species, n_fishery_group] p_observable_sg;

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
  vector[n_species]       deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] deaths_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       cryptic_deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        cryptic_deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] cryptic_deaths_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_method] captures_sm[n_species] = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_method] deaths_sm[n_species] = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_method] cryptic_deaths_sm[n_species] = rep_array(rep_vector(0.0, n_method), n_species);
  
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
  matrix[n_species,n_method]        risk_ratio_sm;
  
  // calculate vulnerability
  {
      real cm;
      real p_live_capture_sg;
      
      for(s in 1:n_species){
        for(g in 1:n_fishery_group){
            
            p_live_capture_sg = p_live_capture_zg[species_group_s[s],g];
            
            // Longline
            if (idx_method_fg[g] == 1 || idx_method_fg[g] == 2) {
                
                cm = cm_longline_par1[s] > 0 ? lognormal_rng(log(cm_longline_par1[s]), cm_longline_par2[s]) : 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg + (1 - p_live_capture_sg) * cm);
                
            } else { 
            
            // Set net
            if (idx_method_fg[g] == 3) {
                
                cm = 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg + (1 - p_live_capture_sg) * cm);
                
            } else { 
            
            // Trawl
            if (idx_method_fg[g] == 4) {
                
                vulnerability_sg[s,g] = q_sg[s,g] * p_live_capture_sg;
                
                // Net
                cm = 1.0 + cm_net_par1[s] > 0 ? lognormal_rng(log(cm_net_par1[s]), cm_net_par2[s]) : 0.0;
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg) * p_capture_type_t[capture_type_group_s[s], 1] * cm;
                
                // Warp
                cm = cm_warp_par1[s] > 0 ? lognormal_rng(log(cm_warp_par1[s]), cm_warp_par2[s]) : 1.0;
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg) * p_capture_type_t[capture_type_group_s[s], 2] * cm;
                
                // Other
                cm = 1.0;
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg) * p_capture_type_t[capture_type_group_s[s], 3] * cm;
                
                // Unknown
                cm = 1.0;
                
                vulnerability_sg[s,g] += q_sg[s,g] * (1 - p_live_capture_sg) * p_capture_type_t[capture_type_group_s[s], 4] * cm;
            }}}
            
            p_observable_sg[s,g] = q_sg[s,g] / vulnerability_sg[s,g];
        }
      }
  }
  
  // Generate priors
  {
      real parA;
      real parB;
      
      // raw values
      for (s in 1:n_species) {
          
        parA = n_breeding_pairs_p1[s];
        parB = n_breeding_pairs_p2[s];

        // N_BP:
        if (n_breeding_pairs_type[s] == 0) {
          prior_n_breeding_pairs_raw_s[s] = uniform_rng(parA, parB);
        }
        if (n_breeding_pairs_type[s] == 1) {
          prior_n_breeding_pairs_raw_s[s] = lognormal_rng(log(parA) - square(parB) / 2, parB);
        }
        if (n_breeding_pairs_type[s] == 2) {
          prior_n_breeding_pairs_raw_s[s] = normal_rng(parA, parB);
        }

        // P_B:
        parA = p_breeding_p1[s];
        parB = p_breeding_p2[s];
        prior_p_breeding_raw_s[s] = normal_rng(logit(parA), parB / (parA * (1.0 - parA)));
      }
      
      // derived values
      for (s in 1:n_species) {

        // N_BP:
        prior_n_breeding_pairs_s[s] = prior_n_breeding_pairs_raw_s[s];

        // P_B: logit-normal
        prior_p_breeding_s[s] = inv_logit(prior_p_breeding_raw_s[s]);

        // N_adults (per annum)
        prior_n_adults_s[s] = 2.0 * prior_n_breeding_pairs_s[s] / prior_p_breeding_s[s];

        // N_vuln_adults (per species and month)
        for(m in 1:n_month){
          prior_n_vuln_adults_s[s, m] = prior_n_adults_s[s] * p_eez[s, m] * (1 - prior_p_breeding_s[s] * p_nest[s, m]);
        }
      }
  }
  
  
  {
      real mu;
      real sigma;
      vector[2] theta;
      vector[1] lambda_max;
        
      for (s in 1:n_species) {
          
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

        // PST
        theta[1] = age_breeding_current_s[s];
        theta[2] = adult_survival_opt_s[s];
        lambda_max = algebra_solver(solve_lambda, y_guess, theta, x_r, x_i);
        r_max_s[s] = log(lambda_max[1]);
        pst_s[s] = 0.5 * psi * r_max_s[s] * n_adults_s[s];
      }
  }

  // get captures and deaths
  // from observed overlap
  for (i in 1:n_i) {
    
    observed_captures_i[i]      = mu_observed_captures_i[i] > 0.0 ? poisson_rng(mu_observed_captures_i[i]) : 0.0;
    observed_live_captures_i[i] = mu_live_captures_i[i] > 0.0 ? poisson_rng(mu_live_captures_i[i]) : 0.0;
    observed_dead_captures_i[i] = mu_dead_captures_i[i] > 0.0 ? poisson_rng(mu_dead_captures_i[i]) : 0.0;
    observed_net_captures_i[i]  = mu_net_captures_i[i] > 0.0 ? poisson_rng(mu_net_captures_i[i]) : 0.0;
    observed_warp_captures_i[i] = mu_warp_captures_i[i] > 0.0 ? poisson_rng(mu_warp_captures_i[i]) : 0.0;
    
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
        
        // simulate observable captures and deaths
        if (mu_observable_captures_j > 0) {
            
            // posterior prediction of average total captures and deaths
            captures_j = mu_observable_captures_j > 0.0 ? poisson_rng(n_years * mu_observable_captures_j) / n_years : 0.0;
            deaths_j   = mu_observable_deaths_j > 0.0 ? poisson_rng(n_years * mu_observable_deaths_j) / n_years : 0.0;
            
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
        mu_deaths_j = mu_captures_j * (1.0 - p_live_captures_j);
        
        // simulate total captures and deaths
        if (mu_captures_j > 0) {
            
            // posterior prediction of average total captures and deaths
            captures_j = mu_captures_j > 0.0 ? poisson_rng(n_years * mu_captures_j) / n_years : 0.0;
            deaths_j   = mu_deaths_j > 0.0 ? poisson_rng(n_years * mu_deaths_j) / n_years : 0.0;
            
            // sum captures
            captures_sm[species_j[j], method_j[j]] += captures_j;
            
            captures_s[species_j[j]]       += captures_j;
            captures_g[fishery_group_j[j]] += captures_j;
            captures_m[method_j[j]]        += captures_j;
            
            // sum deaths
            deaths_sm[species_j[j], method_j[j]] += deaths_j;
            
            deaths_s[species_j[j]]       += deaths_j;
            deaths_g[fishery_group_j[j]] += deaths_j;
            deaths_m[method_j[j]]        += deaths_j;
            
            // expected cryptic deaths 
            cryptic_deaths_sm[species_j[j], method_j[j]] += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            
            cryptic_deaths_s[species_j[j]]       += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            cryptic_deaths_g[fishery_group_j[j]] += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
            cryptic_deaths_m[method_j[j]]        += deaths_j * (1 - p_observable_sg[species_j[j], fishery_group_j[j]]);
        }
      }
  }
  
  // summmary results by species
  for (s in 1:n_species) {
    
    // risk per species
    risk_ratio_s[s]      = deaths_s[s] / pst_s[s];
    
    for (m in 1:n_method) {
      risk_ratio_sm[s,m] = deaths_sm[s,m] / pst_s[s];
    }
  }
  
  // summary results by method
  for (m in 1:n_method) {
    
    // mean risk accross species weighted by species-specific deaths
    risk_ratio_m[m] = sum((to_vector(risk_ratio_sm[,m]) .* deaths_s) / sum(deaths_s));
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
  traces[5] = vector_norm(append_row(beta_live_capture_g, beta_live_capture_z));
  traces[6] = vector_norm(append_row(append_row(n_breeding_pairs_raw_0, n_breeding_pairs_raw_1), n_breeding_pairs_raw_2));
  traces[7] = vector_norm(p_breeding_raw_s);
  
} // end of generated quantities
