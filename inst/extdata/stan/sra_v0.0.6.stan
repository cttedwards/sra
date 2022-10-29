functions {

  // Objective function for allometric (demographic invariant) 
  // solution to lambda
  vector solve_lambda_di(vector y, vector theta, real[] x_r, int[] x_i) {
    
    // age at first breeding
    // (age of maturity + 1)
    real a = theta[1];
    // adult survivorship
    real s = theta[2];
    // allometric constant
    real k = theta[3];
    
    // maximum growth rate
    // (initial value)
    real lambda = y[1];
    
    // objective value
    vector[1] z;
    z[1] = exp(k * pow(a + (s / (lambda - s)), -1)) - lambda;
    return z;
  }
  
  // Objective function for Euler-Lotka 
  // solution to lambda
  vector solve_lambda_el(vector y, vector theta, real[] x_r, int[] x_i) {
    
    // age at first breeding
    real a = theta[1];
    // adult survivorship
    real s = theta[2];
    // survivorship to
    // age at first breeding
    real sa = theta[3];
    // probability of breeding
    real p = theta[4];
    // mean clutch size
    real c = theta[5];
    
    // fecundity
    // (number of females)
    real f = c * 0.5 * p;
    
    // maximum growth rate
    // (initial value)
    real lambda = y[1];
    
    // objective value
    vector[1] z;
    z[1] = s * pow(lambda, a - 1) - pow(lambda, a) + f * sa;
    return z;
  }
  
  // return the Euclidean norm of a vector
  real vector_norm(vector x) {
	    
	    real i = 0.0;
	    
	    for (j in 1:num_elements(x))
	        i += pow(x[j], 2.0);
	        
	    return pow(i, 0.5);
  }

} // end of functions
data {

  // Dimensions
  int n_species;                      // number of species (71)
  int n_month;                        // number of months (12)
  int n_area;                         // number of areas (for output only)
  int n_method;                       // number of fishing methods (4)
  int n_fishery_group;                // number of fishery groups
  int n_species_group;                // number of species groups
  int n_fishery_group_m[n_method];    // number of fishery groups for each method
  int idx_method_fg[n_fishery_group]; // index vector indicating which method applies to which group
  int n_net_group;
  int idx_net_to_species[n_species];
  
  int n_breeding_pairs_type[n_species];
  int n_breeding_pairs_type_0;
  int n_breeding_pairs_type_1;
  int n_breeding_pairs_type_2;
  int adult_survival_opt_type[n_species];
  int adult_survival_opt_type_0;
  int adult_survival_opt_type_1;
  
  // Cryptic mortality multipliers (log-normal)
  vector[n_species] cc_longline_par1;
  vector[n_species] cc_longline_par2;
  vector[n_species] cc_warp_par1;
  vector[n_species] cc_warp_par2;
  vector[n_species] cc_net_par1;
  vector[n_species] cc_net_par2;
  
  // structural groups per species
  int<lower=1,upper=n_species_group> species_group_s[n_species];      // species groups for vulnerability

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
  // S_opt
  vector<lower=0>[n_species] adult_survival_opt_p1;
  vector<lower=0>[n_species] adult_survival_opt_p2;
  // S_egg
  vector<lower=0>[n_species] egg_survival_mult_p1;
  vector<lower=0>[n_species] egg_survival_mult_p2;
  // S_juv
  vector<lower=0>[n_species] juv_survival_mult_p1;
  vector<lower=0>[n_species] juv_survival_mult_p2;
  // clutch size
  vector[n_species] clutch_size_s;
  
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
  
  // Observed captures per event
  int<lower=0> live_captures_i[n_i];
  int<lower=0> dead_captures_i[n_i];
  int<lower=0> net_captures_i[n_i];
  int<lower=0> warp_captures_i[n_i];
  int<lower=0> captures_i[n_i];
  int use_net_captures_i[n_i];
  
  // Observed density overlap
  vector<lower=0>[n_i] overlap_i;
  
  // All fishing events
  int<lower=1>                       n_j; 
  int<lower=1,upper=n_month>         month_j[n_j];
  int<lower=1,upper=n_area>          area_j[n_j];
  int<lower=1,upper=n_method>        method_j[n_j];
  int<lower=1,upper=n_fishery_group> fishery_group_j[n_j];
  int<lower=1,upper=n_species>       species_j[n_j];
  int<lower=1,upper=n_species_group> species_group_j[n_j];
  //int use_net_captures_j[n_j];
  
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
  vector[1] y_guess = rep_vector(1, 1);
  //y_guess[1] = 1.0;
  
  // Lambda estimation flag
  // 1: DIM
  // 2: Euler-Lotka
  // 3: Dillingham method
  int lambda_method = 3;

} // end of transformed data
parameters {

  // catchability predictors
  vector[n_method]                   log_q_intercept;
  vector[n_fishery_group - n_method] log_q_method_raw;
  vector[n_species_group - 1]        log_q_species_raw;
  
  // catchability
  // random effect
  real<lower=0> tau;
  matrix[n_fishery_group,n_species] eps_gs;

  // Predictors of live capture
  // probability
  vector[n_fishery_group] beta_live_capture_g;
  vector[n_species_group] beta_live_capture_z;
  
  // Predictors of biological values
  vector<lower=0,upper=1>[n_breeding_pairs_type_0] n_breeding_pairs_raw_0; 
  vector<lower=0>[n_breeding_pairs_type_1] n_breeding_pairs_raw_1;
  vector<lower=0>[n_breeding_pairs_type_2] n_breeding_pairs_raw_2;
  vector[n_species] p_breeding_raw_s;
  vector<lower=0,upper=1>[adult_survival_opt_type_0] adult_survival_opt_raw_0; 
  vector<lower=0>[adult_survival_opt_type_1] adult_survival_opt_raw_1;
  vector<lower=0,upper=1>[n_species] age_breeding_current_raw_s;
  
  vector<lower=0>[n_species] egg_survival_mult_raw_s;
  vector<lower=0>[n_species] juv_survival_mult_raw_s;
  
  // Distribution of trawl captures across capture types
  // per capture type group
  vector<lower=0,upper=1>[n_net_group] p_net_capture_z;
  
  vector<lower=1>[n_species] lambda_max[3];
  
} // end of parameters
transformed parameters {

  // catchability (probability of
  // observing a capture)
  vector<lower=0>[n_fishery_group] q_sg[n_species];
  
  // catchability predictands
  vector[n_species_group] log_q_z;
  vector[n_fishery_group] log_q_g;
  
  vector<lower=0,upper=1>[n_species] p_net_capture_s;
  
  // probabilities of live capture
  vector<lower=0,upper=1>[n_fishery_group] p_live_capture_zg[n_species_group];
  vector<lower=0,upper=1>[n_i] p_live_capture_i;
  vector<lower=0,upper=1>[n_i] p_net_capture_i;

  // biological parameters
  vector<lower=0>[n_species] n_breeding_pairs_s;
  vector<lower=0>[n_species] n_adults_s;
  vector<lower=0,upper=1>[n_species] p_breeding_s;
  vector[n_species] age_breeding_current_s;
  vector[n_species] adult_survival_opt_s;
  vector[n_species] egg_survival_opt_s;
  vector[n_species] juv_survival_opt_s;
  vector[n_species] juv_survival_to_adult_opt_s;
  vector<lower=0>[n_month]   n_vuln_adults_s[n_species];
  
  // expected observed captures
  // for diagnosing model fit
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
    
    // net capture by species
    for (s in 1:n_species) {
        p_net_capture_s[s] = p_net_capture_z[idx_net_to_species[s]];
    }
    
    // back transform catchability onto natural scale
    for (g in 1:n_fishery_group) {
      for (s in 1:n_species) {
        q_sg[s,g] = exp(log_q_z[species_group_s[s]] + log_q_g[g] + eps_gs[g,s]);
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
      
      int c0 = 1;
      int c1 = 1;
      
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
        
        // S_opt - uniform or logit-normal
        parA = adult_survival_opt_p1[s];
        parB = adult_survival_opt_p2[s];
        
        if (adult_survival_opt_type[s] == 0) {
          adult_survival_opt_s[s] = parA + (parB - parA) * adult_survival_opt_raw_0[c0];
          c0 += 1;
        } else {
          adult_survival_opt_s[s] = inv_logit(adult_survival_opt_raw_1[c1]);
          c1 += 1;
        }
        
        // S_egg
        egg_survival_opt_s[s] = adult_survival_opt_s[s] * inv_logit(egg_survival_mult_raw_s[s]);
        
        // S_juv
        juv_survival_opt_s[s] = adult_survival_opt_s[s] * inv_logit(juv_survival_mult_raw_s[s]);
        
        // A_curr: current age at first breeding 
        parA = age_breeding_current_p1[s];
        parB = age_breeding_current_p2[s];
        age_breeding_current_s[s] = parA + (parB - parA) * age_breeding_current_raw_s[s];
        
        // S_to_adult
        juv_survival_to_adult_opt_s[s] = egg_survival_opt_s[s] * pow(juv_survival_opt_s[s], age_breeding_current_s[s] - 1.0);
        
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
          
        // density overlap
        density_overlap_i = fmax(1e-8, overlap_i[i] * n_vuln_adults_s[species_i[i], month_i[i]]);
        
        // expected observed captures
        mu_observed_captures_i[i] = q_sg[species_i[i],fishery_group_i[i]] * density_overlap_i;
        
        // probability that a trawl capture is a net capture
        p_net_capture_i[i] = use_net_captures_i[i] > 0 ? p_net_capture_s[species_i[i]] : 1.0;
        
        // conditional probability of live capture
        // (conditional on it being a net capture
        // for the trawl fishery)
        //p_conditional_live_capture_i[i] = p_live_capture_zg[species_group_i[i], fishery_group_i[i]];
        
        // conditional probability of live capture
        p_live_capture_i[i] = p_live_capture_zg[species_group_i[i], fishery_group_i[i]];
        
        // expected live and dead captures
        // (for trawl this accounts for the fact
        // that only net captures can be alive)
        mu_live_captures_i[i] = mu_observed_captures_i[i] * p_live_capture_i[i];
        mu_dead_captures_i[i] = mu_observed_captures_i[i] * (1.0 - p_live_capture_i[i]);

        // expected net and warp captures
        // (scaled to match number of observations)
        mu_net_captures_i[i]  = use_net_captures_i[i] > 0 ? mu_observed_captures_i[i] * p_net_capture_i[i] * (net_captures_i[i] + warp_captures_i[i]) / captures_i[i] : 0.0;
        mu_warp_captures_i[i] = use_net_captures_i[i] > 0 ? mu_observed_captures_i[i] * (1.0 - p_net_capture_i[i]) * (net_captures_i[i] + warp_captures_i[i]) / captures_i[i] : 0.0;
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
  to_vector(eps_gs) ~ normal(0.0, tau);
  tau ~ cauchy(0, 1);
  
  // lambda is an estimatd parameter
  {
      vector[5] theta = rep_vector(0.0, 5);
      real prior_lambda_max;
        
      for (s in 1:n_species) {

        // Lambda (Demographic Invariant)
        theta[1] = age_breeding_current_s[s];
        theta[2] = adult_survival_opt_s[s];
        theta[3] = 1.0;
        prior_lambda_max = algebra_solver(solve_lambda_di, y_guess, theta, x_r, x_i)[1];
        
        lambda_max[1,s] ~ normal(prior_lambda_max, 0.01);
        lambda_max[3,s] ~ normal(prior_lambda_max, 0.01);
        
        // Lambda (Euler-Lotka)
        theta[1] = age_breeding_current_s[s];
        theta[2] = adult_survival_opt_s[s];
        theta[3] = juv_survival_to_adult_opt_s[s];
        theta[4] = p_breeding_s[s];
        theta[5] = clutch_size_s[s];
        prior_lambda_max = algebra_solver(solve_lambda_el, y_guess, theta, x_r, x_i)[1];
        
        lambda_max[2,s] ~ normal(prior_lambda_max, 0.01);
        lambda_max[3,s] ~ normal(prior_lambda_max, 0.01);
      }
  }
  
  // priors on estimated biological 
  // parameters
  {
      real parA;
      real parB;
      
      int k0 = 1;
      int k1 = 1;
      int k2 = 1;
      
      int c0 = 1;
      int c1 = 1;
      
      for (s in 1:n_species) {
          
        // N_BP:
        parA = n_breeding_pairs_p1[s];
        parB = n_breeding_pairs_p2[s];
        
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
        
        // S_opt:
        parA = adult_survival_opt_p1[s];
        parB = adult_survival_opt_p2[s];
        
        if (adult_survival_opt_type[s] == 0) {
          adult_survival_opt_raw_0[c0] ~ uniform(0, 1);
          c0 += 1;
        }
        if (adult_survival_opt_type[s] == 1) {
          adult_survival_opt_raw_1[c1] ~ normal(logit(parA), parB / (parA * (1.0 - parA)));
          c1 += 1;
        }
        
        // S_egg:
        parA = egg_survival_mult_p1[s];
        parB = egg_survival_mult_p2[s];
        egg_survival_mult_raw_s[s] ~ normal(logit(parA), parB / (parA * (1.0 - parA)));
        
        // S_juv:
        parA = juv_survival_mult_p1[s];
        parB = juv_survival_mult_p2[s];
        juv_survival_mult_raw_s[s] ~ normal(logit(parA), parB / (parA * (1.0 - parA)));
        
        // A_curr
        age_breeding_current_raw_s ~ uniform(0, 1);
        
        // P_B:
        parA = p_breeding_p1[s];
        parB = p_breeding_p2[s];
        p_breeding_raw_s[s] ~ normal(logit(parA), parB / (parA * (1.0 - parA)));
      }
  }

  // priors on live capture coefficients
  beta_live_capture_g ~ std_normal();
  beta_live_capture_z ~ std_normal();

  // likelihood for total observed captures
  captures_i ~ poisson(mu_observed_captures_i);
  
  // likelihood for net capture
  // given trawl capture with recorded
  // trawl type
  for (i in 1:n_i) {
      if (use_net_captures_i[i] > 0) {
          net_captures_i[i] ~ binomial(net_captures_i[i] + warp_captures_i[i], p_net_capture_i[i]);
      }
  }
  
  // likelihood for live capture
  // given capture (for the trawl
  // fishery all live captures are
  // assumed to be net captures)
  for (i in 1:n_i) {     
      if (captures_i[i] > 0) {
          live_captures_i[i] ~ binomial(captures_i[i], p_live_capture_i[i]);
      }
  }
} // end of model
generated quantities {

  // summary traces
  vector[11] traces;
  
  // Priors for biological parameters
  vector[n_species] prior_n_breeding_pairs_raw_s; 
  vector[n_species] prior_p_breeding_raw_s;
  
  // priors for derived biological parameters
  vector<lower=0>[n_species] prior_n_breeding_pairs_s;
  vector<lower=0>[n_species] prior_n_adults_s;
  vector<lower=0>[n_month]   prior_n_vuln_adults_s[n_species];
  vector<lower=0,upper=1>[n_species] prior_p_breeding_s;

  // posterior predictive checking
  int observed_captures_i[n_i];
  int observed_live_captures_i[n_i];
  int observed_dead_captures_i[n_i];
  int observed_net_captures_i[n_i];
  int observed_warp_captures_i[n_i];
  
  // back calculation of vulnerability
  // (probability of capture)
  matrix[n_species, n_fishery_group] vulnerability_sg;
  matrix[n_species_group, n_fishery_group] vulnerability_zg = rep_matrix(0.0, n_species_group, n_fishery_group);
  matrix[n_species_group, n_method]        vulnerability_zm = rep_matrix(0.0, n_species_group, n_method);
  // (catchability)
  matrix[n_species_group, n_fishery_group] q_zg = rep_matrix(0.0, n_species_group, n_fishery_group);
  matrix[n_species_group, n_method]        q_zm = rep_matrix(0.0, n_species_group, n_method);
  
  // record overal cryptic capture multipliers
  matrix[n_species, n_fishery_group]       cryptic_multiplier_sg;
  matrix[n_species_group, n_fishery_group] cryptic_multiplier_zg = rep_matrix(0.0, n_species_group, n_fishery_group);
  matrix[n_species_group, n_method]        cryptic_multiplier_zm = rep_matrix(0.0, n_species_group, n_method);

  // Captures and deaths
  vector[n_species]       observed_captures_s = rep_vector(0.0, n_species);
  vector[n_method]        observed_captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observed_captures_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       observed_net_captures_s = rep_vector(0.0, n_species);
  vector[n_method]        observed_net_captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observed_net_captures_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       observed_dead_captures_s = rep_vector(0.0, n_species);
  vector[n_method]        observed_dead_captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] observed_dead_captures_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_species]       captures_s = rep_vector(0.0, n_species);
  vector[n_method]        captures_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] captures_g = rep_vector(0.0, n_fishery_group);
  
  vector[n_species]       deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] deaths_g = rep_vector(0.0, n_fishery_group);
  vector[n_species]       cryptic_deaths_s = rep_vector(0.0, n_species);
  vector[n_method]        cryptic_deaths_m = rep_vector(0.0, n_method);
  vector[n_fishery_group] cryptic_deaths_g = rep_vector(0.0, n_fishery_group);
  
  real p_survive = uniform_rng(0, 1);
  
  vector[n_method]        captures_sm[n_species]            = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_fishery_group] captures_sg[n_species]            = rep_array(rep_vector(0.0, n_fishery_group), n_species);
  vector[n_fishery_group] captures_zg[n_species_group]      = rep_array(rep_vector(0.0, n_fishery_group), n_species_group);
  vector[n_method]        dead_captures_sm[n_species]       = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_fishery_group] dead_captures_zg[n_species_group] = rep_array(rep_vector(0.0, n_fishery_group), n_species_group);
  
  vector[n_method]        deaths_sm[n_species]                = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_fishery_group] deaths_sg[n_species]                = rep_array(rep_vector(0.0, n_fishery_group), n_species);
  vector[n_fishery_group] deaths_zg[n_species_group]          = rep_array(rep_vector(0.0, n_fishery_group), n_species_group);
  vector[n_method]        cryptic_deaths_sm[n_species]        = rep_array(rep_vector(0.0, n_method), n_species);
  vector[n_fishery_group] cryptic_deaths_sg[n_species]        = rep_array(rep_vector(0.0, n_fishery_group), n_species);
  vector[n_fishery_group] cryptic_deaths_zg[n_species_group]  = rep_array(rep_vector(0.0, n_fishery_group), n_species_group);
  
  // rmax, PST, risk ratio
  vector[n_species] r_max_s;
  vector[n_species] pst_s;
  
  vector[n_species]          risk_ratio_s;
  vector[n_species]          pop_status_s;
  vector[n_method]           risk_ratio_m;
  matrix[n_species,n_method] risk_ratio_sm;
  
  // calculate vulnerability
  {
      real cc_ll;
      real cc_sn;
      real cc_net;
      real cc_warp;
      
      real omega = uniform_rng(0, 1);
      
      real p_live_capture_sg;
      
      for(s in 1:n_species){
        for(g in 1:n_fishery_group){
            
            // probability of live capture
            // (conditional live capture for trawl fishery)
            p_live_capture_sg = p_live_capture_zg[species_group_s[s],g];
            
            // Longline
            // (dead captures are cryptic)
            if (idx_method_fg[g] == 1 || idx_method_fg[g] == 2) {
                
                cc_ll = cc_longline_par1[s] > 0.0 ? lognormal_rng(log(cc_longline_par1[s]) - 0.5 * pow(cc_longline_par2[s], 2), cc_longline_par2[s]) : 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg + (1 - p_live_capture_sg) * cc_ll);
                
                cryptic_multiplier_sg[s,g] = p_live_capture_sg * (1 - omega) + (1 - p_live_capture_sg) * cc_ll;
                
            } else { 
            
            // Set net
            // (no captures are cryptic)
            if (idx_method_fg[g] == 3) {
                
                cc_sn = 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_live_capture_sg + (1 - p_live_capture_sg) * cc_sn);
                
                cryptic_multiplier_sg[s,g] = p_live_capture_sg * (1 - omega) + (1 - p_live_capture_sg) * cc_sn;
                
            } else { 
            
            // Trawl
            if (idx_method_fg[g] == 4) {
                
                // Net
                cc_net = cc_net_par1[s] > 0.0 ? 1.0 + lognormal_rng(log(cc_net_par1[s]) - 0.5 * pow(cc_net_par2[s], 2), cc_net_par2[s]) : 0.0;
                
                // Warp
                cc_warp = cc_warp_par1[s] > 0.0 ? lognormal_rng(log(cc_warp_par1[s]) - 0.5 * pow(cc_warp_par2[s], 2), cc_warp_par2[s]) : 1.0;
                
                vulnerability_sg[s,g] = q_sg[s,g] * (p_net_capture_s[s] * cc_net + (1.0 - p_net_capture_s[s]) * cc_warp);
                
                cryptic_multiplier_sg[s,g] = p_net_capture_s[s] * cc_net * (1 - p_live_capture_sg * omega) + (1.0 - p_net_capture_s[s]) * cc_warp;
                
            }}}
        }
      }
      
      // Assign geometric or arithmetic mean to species group
      {
          int n[n_species_group] = rep_array(0, n_species_group);
          
          for(s in 1:n_species){
              for(g in 1:n_fishery_group){
                  
                q_zg[species_group_s[s],g]                  += log(q_sg[s,g]);
                vulnerability_zg[species_group_s[s],g]      += log(vulnerability_sg[s,g]);

                cryptic_multiplier_zg[species_group_s[s],g] += cryptic_multiplier_sg[s,g];
                
              }
              
              // number of species per group
              n[species_group_s[s]] += 1;
          }
          
          for(z in 1:n_species_group){
              for(g in 1:n_fishery_group){
                  
                q_zg[z,g]                   = exp(inv(n[z]) * q_zg[z,g]);
                vulnerability_zg[z,g]       = exp(inv(n[z]) * vulnerability_zg[z,g]);

                cryptic_multiplier_zg[z,g] /= n[z];
              }
          }
      }
      
      // Assign mean to species group and method
      {
          int n[n_species_group, n_method] = rep_array(0, n_species_group, n_method);
          
          for(s in 1:n_species){
              for(g in 1:n_fishery_group){
                  
                q_zm[species_group_s[s],idx_method_fg[g]]                  += log(q_sg[s,g]);
                vulnerability_zm[species_group_s[s],idx_method_fg[g]]      += log(vulnerability_sg[s,g]);
                
                cryptic_multiplier_zm[species_group_s[s],idx_method_fg[g]] += cryptic_multiplier_sg[s,g];
                
                n[species_group_s[s],idx_method_fg[g]] += 1;
              }
          }
          
          for(z in 1:n_species_group){
              for(m in 1:n_method){
                  
                q_zm[z,m]                   = exp(inv(n[z,m]) * q_zm[z,m]);
                vulnerability_zm[z,m]       = exp(inv(n[z,m]) * vulnerability_zm[z,m]);
                
                cryptic_multiplier_zm[z,m] /= n[z,m];
              }
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
  
  for (s in 1:n_species) {

     // calculate PST
     r_max_s[s] = log(lambda_max[lambda_method,s]);
     pst_s[s] = 0.5 * psi * r_max_s[s] * n_adults_s[s];
  }

  // get posterior prediction of
  // observed captures
  // from observed overlap
  for (i in 1:n_i) {
  
    observed_captures_i[i]          = 0;
    observed_live_captures_i[i]     = 0;
    observed_dead_captures_i[i]     = 0;
    observed_net_captures_i[i]      = 0;
    observed_warp_captures_i[i]     = 0;
    
    if (mu_observed_captures_i[i] > 0) {
    
        observed_captures_i[i] = poisson_rng(mu_observed_captures_i[i]);
        
        observed_live_captures_i[i] = poisson_rng(mu_live_captures_i[i]);
        observed_dead_captures_i[i] = poisson_rng(mu_dead_captures_i[i]);
       
        if (use_net_captures_i[i] > 0) {
            
            observed_net_captures_i[i]  = poisson_rng(mu_net_captures_i[i]);
            observed_warp_captures_i[i] = poisson_rng(mu_warp_captures_i[i]);
        }
    }    
        
    // sum captures
    observed_captures_s[species_i[i]]       += observed_captures_i[i];
    observed_captures_g[fishery_group_i[i]] += observed_captures_i[i];
    observed_captures_m[method_i[i]]        += observed_captures_i[i];
    
    // sum net captures
    observed_net_captures_s[species_i[i]]       += observed_net_captures_i[i];
    observed_net_captures_g[fishery_group_i[i]] += observed_net_captures_i[i];
    observed_net_captures_m[method_i[i]]        += observed_net_captures_i[i];
    
    // sum dead captures
    observed_dead_captures_s[species_i[i]]       += observed_dead_captures_i[i];
    observed_dead_captures_g[fishery_group_i[i]] += observed_dead_captures_i[i];
    observed_dead_captures_m[method_i[i]]        += observed_dead_captures_i[i];
  }

  // get annual captures and deaths
  // from total overlap
  {
    real density_overlap_j;
    real captures_j;
    real net_captures_j;
    real deaths_j;
    real observable_deaths_j;
    
    for (j in 1:n_j) {

        density_overlap_j = fmax(1e-8, overlap_j[j] * n_vuln_adults_s[species_j[j], month_j[j]]);

        // observable captures per year
        captures_j = poisson_rng(q_sg[species_j[j], fishery_group_j[j]] * density_overlap_j) / n_years;

        // death
        deaths_j = captures_j * cryptic_multiplier_sg[species_j[j], fishery_group_j[j]];
        
        // observable death
        observable_deaths_j = captures_j * (1 - p_live_capture_zg[species_group_j[j], fishery_group_j[j]]);
        
        // sum captures
        captures_sm[species_j[j], method_j[j]]                         += captures_j;
        captures_sg[species_j[j], fishery_group_j[j]]                  += captures_j;
        captures_zg[species_group_s[species_j[j]], fishery_group_j[j]] += captures_j;
            
        captures_s[species_j[j]]       += captures_j;
        captures_g[fishery_group_j[j]] += captures_j;
        captures_m[method_j[j]]        += captures_j;

        // sum deaths
        deaths_sm[species_j[j], method_j[j]]                         += deaths_j;
        deaths_sg[species_j[j], fishery_group_j[j]]                  += deaths_j;
        deaths_zg[species_group_s[species_j[j]], fishery_group_j[j]] += deaths_j;

        deaths_s[species_j[j]]       += deaths_j;
        deaths_g[fishery_group_j[j]] += deaths_j;
        deaths_m[method_j[j]]        += deaths_j;

        // expected cryptic deaths 
        cryptic_deaths_sm[species_j[j], method_j[j]]                         += deaths_j - observable_deaths_j;
        cryptic_deaths_sg[species_j[j], fishery_group_j[j]]                  += deaths_j - observable_deaths_j;
        cryptic_deaths_zg[species_group_s[species_j[j]], fishery_group_j[j]] += deaths_j - observable_deaths_j;

        cryptic_deaths_s[species_j[j]]       += deaths_j - observable_deaths_j;
        cryptic_deaths_g[fishery_group_j[j]] += deaths_j - observable_deaths_j;
        cryptic_deaths_m[method_j[j]]        += deaths_j - observable_deaths_j;
   }
  }
  
  // summmary results by species
  for (s in 1:n_species) {
    
    // risk per species
    risk_ratio_s[s] = deaths_s[s] / pst_s[s];
    pop_status_s[s] = 1 - deaths_s[s] / (n_adults_s[s] * r_max_s[s]);
    
    for (m in 1:n_method) {
      risk_ratio_sm[s,m] = deaths_sm[s,m] / pst_s[s];
    }
  }
  
  // summary results by method
  for (m in 1:n_method) {
    
    // mean risk across species weighted by species-specific deaths
    risk_ratio_m[m] = sum((to_vector(risk_ratio_sm[,m]) .* deaths_s) / sum(deaths_s));
  }
  
  // Mortality constraint
  // **
  
  // summary trace outputs
  traces[1] = vector_norm(log_q_intercept);
  traces[2] = vector_norm(log_q_method_raw);
  traces[3] = vector_norm(log_q_species_raw);
  traces[4] = tau;
  traces[5] = vector_norm(append_row(beta_live_capture_g, beta_live_capture_z));
  traces[6] = vector_norm(append_row(append_row(n_breeding_pairs_raw_0, n_breeding_pairs_raw_1), n_breeding_pairs_raw_2));
  traces[7] = vector_norm(p_breeding_raw_s);
  traces[8] = vector_norm(append_row(adult_survival_opt_raw_0, adult_survival_opt_raw_1));
  traces[9] = vector_norm(age_breeding_current_raw_s);
  traces[10] = vector_norm(append_row(egg_survival_mult_raw_s, juv_survival_mult_raw_s));
  traces[11] = vector_norm(p_net_capture_z);
  
} // end of generated quantities
