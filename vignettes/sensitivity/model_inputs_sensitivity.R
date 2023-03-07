
##################################
# MANUAL RUN OF DATA PREPARATION #
# CODE FOR MODEL SENSITIVITIES   #
# +> ACCESS DATA                 #
# DIRECTLY                       #
##################################

# DATA DIRECTORY
dat_dir <- file.path("..", paste0("data-v", packageVersion("sraInputsBio")))
dir.create(dat_dir)

# INPUT DIR
input_dir <- file.path("../../sraInputs", paste0("data-v", packageVersion("sraInputsBio")))

library(dplyr)
library(ggplot2)

# Get Version
VERSION <- package_version(scan('../DESCRIPTION',what = character(),skip = 3,nlines = 1)[2])

# from manual run of sraInputs package data prep
load(file.path(input_dir, "methods.rda"))
load(file.path(input_dir, "fishing_groups.rda"))
load(file.path(input_dir, "fishing_years_risk.rda"))
load(file.path(input_dir, "overlap_o.rda"))
load(file.path(input_dir, "overlap_t.rda"))

library(sraUtils)

# from sraInputs package build
data("species", package = "sraInputs")
data("species_groups", package = "sraInputs")
data("cc_net", package = "sraInputs")
data("cc_warp", package = "sraInputs")
data("cc_longline", package = "sraInputs")
data("fishing_group_definitions", package = "sraInputs")

library(sraInputsBio)

# get fecundity parameters
data("inputsBio_v2.0.0", package = "sraInputsBio")

stopifnot(isTRUE(all.equal(attributes(inputsBio_v2.0.0)$version, packageVersion("sraInputsBio"))))

S_egg  <- inputsBio_v2.0.0[["S_egg_mult"]]
S_juv  <- inputsBio_v2.0.0[["S_juv_mult"]]
clutch <- inputsBio_v2.0.0[["clutch_size"]][,"clutch_size"]

# other data
data_file <- paste0("inputsBio_v", packageVersion("sraInputsBio"))
do.call(data, list(data_file, package = "sraInputsBio"))

# get data
inputsBio <- do.call(get, list(data_file))
rm(data_file)

# check that package data and version match
stopifnot(isTRUE(all.equal(attributes(inputsBio)$version, packageVersion("sraInputsBio"))))

# Import demographic data
N_BP   <- inputsBio[["N_BP"]]
S_curr <- inputsBio[["S_curr"]]
S_opt  <- inputsBio[["S_opt"]]
A_curr <- inputsBio[["A_curr"]]
P_B    <- inputsBio[["P_B"]]
p_nest <- inputsBio[["p_nest"]]
p_eez  <- inputsBio[["p_EEZ"]]

p_survive <- matrix(1.0, nrow = length(species), ncol = length(unique(methods$id_method)))

# log-uniform parameters should be uniform
N_BP$distribution[grepl("log-uniform", N_BP$distribution)] <- "uniform"

N_BP   <- N_BP   %>% filter(code_sra %in% species)
S_curr <- S_curr %>% filter(code_sra %in% species)
S_opt  <- S_opt  %>% filter(code_sra %in% species)
S_egg  <- S_egg  %>% filter(code_sra %in% species)
S_juv  <- S_juv  %>% filter(code_sra %in% species)
A_curr <- A_curr %>% filter(code_sra %in% species)
P_B    <- P_B    %>% filter(code_sra %in% species)

p_eez  <- p_eez  %>% subset(rownames(.) %in% species)
p_nest <- p_nest %>% subset(rownames(.) %in% species)


# Check that there are always some birds in zone when there are non-zero captures
tab1 <- tapply(overlap_o$n_captures, list(overlap_o$species, overlap_o$month), sum, na.rm = TRUE)
tab1 <- tab1[match(species_groups$code_sra, rownames(tab1)),]

tab2 <- tapply(overlap_o$overlap, list(overlap_o$species, overlap_o$month), sum, na.rm = TRUE)
tab2 <- tab2[match(species_groups$code_sra, rownames(tab2)),]

tab <- full_join(reshape2::melt(tab1, varnames = c("species", "month")), reshape2::melt(tab2, varnames = c("species", "month")), by = c("species", "month"))
stopifnot(nrow(tab) == length(species) * 12)
tab <- tab %>% filter(!is.na(value.x))

stopifnot(with(tab, !any(value.x > 0 & value.y == 0)))


# check that p_nest > 0 for breeding season
#stopifnot(all(p_nest[breeding > 0] > 0))

# check that p_nest == 0 for non-breeding season
#stopifnot(all(p_nest[breeding == 0] == 0))

# check data are orthogonal
stopifnot(all(with(overlap_o, table(species, id_method, month, id_fishing_group)) %in% 0:1))
stopifnot(all(with(overlap_t, table(species, id_method, month, id_fishing_group)) %in% 0:1))


# MODEL DATA
sra_dat <- list(n_species = length(species),
            n_month = 12L,
            n_area = 1L,
            n_method = nrow(methods),
            n_species_group = length(unique(species_groups$species_group)),
            n_fishery_group = length(unique(fishing_groups$id_fishing_group)), 
            n_fishery_group_m = as.integer(table(fishing_groups$method)),
            n_net_group = length(unique(cc_warp$par1)),
            
            species_group_s = species_groups$id_species_group,
            
            idx_method_fg = fishing_groups$id_method,
            idx_net_to_species = match(cc_warp$par1, unique(cc_warp$par1)),
            
            cc_longline_par1 = cc_longline$par1,
            cc_longline_par2 = cc_longline$par2,
            cc_warp_par1 = cc_warp$par1,
            cc_warp_par2 = cc_warp$par2,
            cc_net_par1 = cc_net$par1,
            cc_net_par2 = cc_net$par2,

            n_breeding_pairs_type = ifelse(N_BP$distribution == 'uniform', 0, ifelse(N_BP$distribution == 'log-normal', 1, 2)),
            n_breeding_pairs_p1 = N_BP$par1,
            n_breeding_pairs_p2 = N_BP$par2,
            n_breeding_pairs_type_0 = sum(N_BP$distribution == 'uniform'),
            n_breeding_pairs_type_1 = sum(N_BP$distribution == 'log-normal'),
            n_breeding_pairs_type_2 = sum(N_BP$distribution == 'normal'),
            p_breeding_p1 = P_B$par1,
            p_breeding_p2 = P_B$par2,
            age_breeding_current_p1 = A_curr$par1,
            age_breeding_current_p2 = A_curr$par2,
            adult_survival_current_type = ifelse(S_curr$distribution == 'uniform', 0, 1),
            adult_survival_current_p1 = S_curr$par1,
            adult_survival_current_p2 = S_curr$par2,
            adult_survival_opt_type = ifelse(S_opt$distribution == 'uniform', 0, 1),
            adult_survival_opt_type_0 = sum(S_opt$distribution == 'uniform'),
            adult_survival_opt_type_1 = sum(S_opt$distribution == 'logit-normal'),
            adult_survival_opt_p1 = S_opt$par1,
            adult_survival_opt_p2 = S_opt$par2,
            egg_survival_mult_p1 = S_egg$par1,
            egg_survival_mult_p2 = S_egg$par2,
            juv_survival_mult_p1 = S_juv$par1,
            juv_survival_mult_p2 = S_juv$par2,
            
            clutch_size_s = clutch,
            
            p_eez  = p_eez,
            p_nest = p_nest,
             
            n_i              = nrow(overlap_o),
            month_i          = overlap_o$month,
            method_i         = overlap_o$id_method,
            fishery_group_i  = overlap_o$id_fishing_group,
            species_i        = overlap_o$id_species,
            species_group_i  = overlap_o$id_species_group,
            overlap_i        = overlap_o$overlap,
            live_captures_i  = overlap_o$n_captures_alive,
            dead_captures_i  = overlap_o$n_captures_dead,
            net_captures_i   = overlap_o$n_captures_net,
            warp_captures_i  = overlap_o$n_captures_warp,
            captures_i       = overlap_o$n_captures,
            use_net_captures_i = overlap_o$use_captures_net,
 
            n_j             = nrow(overlap_t),
            month_j         = overlap_t$month,
            area_j          = rep(1L, nrow(overlap_t)),
            method_j        = overlap_t$id_method,
            fishery_group_j = overlap_t$id_fishing_group,
            species_j       = overlap_t$id_species,
            species_group_j = overlap_t$id_species_group,
            overlap_j       = overlap_t$overlap,
            
            n_years = length(fishing_years_risk),
            
            psi = 1.0,
            
            lambda_method = 3L)


attributes(sra_dat)$version <- VERSION


# check all probability distributions are represented in the model code
stopifnot(isTRUE(all(sra_dat$n_breeding_pairs_type %in% 0:2)))
stopifnot(isTRUE(all(sra_dat$adult_survival_current_type %in% 0:1)))
stopifnot(isTRUE(all(sra_dat$adult_survival_opt_type %in% 0:1)))


# initial values
sra_ini <- list(log_q_intercept = rep(-10, sra_dat$n_method),
                log_q_method_raw = rep(0, sra_dat$n_fishery_group - sra_dat$n_method),
                log_q_species_raw = rep(0.0, sra_dat$n_species_group - 1),
                tau = 1.0,
                eps_gs = matrix(0.0, sra_dat$n_fishery_group, sra_dat$n_species),
                beta_live_capture_fg = rep(2.0, sra_dat$n_fishery_group),
                beta_live_capture_sg = rep(0.0, sra_dat$n_species_group),
                n_breeding_pairs_raw_0 = structure(rep(0.5, sra_dat$n_breeding_pairs_type_0), dim = sra_dat$n_breeding_pairs_type_0),
                n_breeding_pairs_raw_1 = structure(with(subset(N_BP, distribution == 'log-normal'), par1), dim = sra_dat$n_breeding_pairs_type_1),
                n_breeding_pairs_raw_2 = structure(with(subset(N_BP, distribution == 'normal'), par1), dim = sra_dat$n_breeding_pairs_type_2),
                p_breeding_raw_s = rep(2.0, sra_dat$n_species),
                adult_survival_opt_raw_0 = rep(0.5, sra_dat$adult_survival_opt_type_0),
                adult_survival_opt_raw_1 = logit(sra_dat$adult_survival_opt_p1[sra_dat$adult_survival_opt_type == 1]),
                egg_survival_mult_raw_s = rep(0.7, sra_dat$n_species),
                juv_survival_mult_raw_s = rep(0.7, sra_dat$n_species),
                age_breeding_current_raw_s = rep(0.5, sra_dat$n_species),
                p_net_capture_z = rep(0.5, sra_dat$n_net_group))
                #p_net_capture_z = rep(0.0, sra_dat$n_species_group),
                #p_net_capture_g = rep(0.0, sum(sra_dat$idx_trawl_fg)))


attributes(sra_ini)$version <- VERSION


# save
save(methods,        file = file.path(dat_dir, "methods.rda"))
save(fishing_groups, file = file.path(dat_dir, "fishing_groups.rda"))
save(sra_dat,        file = file.path(dat_dir, "sra_dat.rda"))


# Building model ----

sra_mdl <- rstan::stan_model(file = "../inst/extdata/stan/sra_v0.0.6.stan")

#model <- sra::sra()

# Get initial values ----

tmp <- rstan::optimizing(object = sra_mdl, data = sra_dat, init = sra_ini, verbose = TRUE, as_vector = FALSE)

sra_ini$log_q_intercept           <- tmp$par$log_q_intercept
sra_ini$log_q_method_raw          <- tmp$par$log_q_method_raw
sra_ini$log_q_species_raw         <- tmp$par$log_q_species_raw
sra_ini$eps_gs                    <- tmp$par$eps_gs
sra_ini$beta_live_capture_g       <- tmp$par$beta_live_capture_g
sra_ini$beta_live_capture_z       <- tmp$par$beta_live_capture_z
sra_ini$n_breeding_pairs_raw_0    <- structure(tmp$par$n_breeding_pairs_raw_0, dim = sra_dat$n_breeding_pairs_type_0)
sra_ini$n_breeding_pairs_raw_1    <- structure(tmp$par$n_breeding_pairs_raw_1, dim = sra_dat$n_breeding_pairs_type_1)
sra_ini$n_breeding_pairs_raw_2    <- structure(tmp$par$n_breeding_pairs_raw_2, dim = sra_dat$n_breeding_pairs_type_2)
sra_ini$p_breeding_raw_s          <- tmp$par$p_breeding_raw_s
sra_ini$adult_survival_opt_raw_0  <- structure(tmp$par$adult_survival_opt_raw_0, dim = sra_dat$adult_survival_opt_type_0)
sra_ini$adult_survival_opt_raw_1  <- structure(tmp$par$adult_survival_opt_raw_1, dim = sra_dat$adult_survival_opt_type_1)
sra_ini$egg_survival_mult_raw_s   <- tmp$par$egg_survival_mult_raw_s
sra_ini$juv_survival_mult_raw_s   <- tmp$par$juv_survival_mult_raw_s
sra_ini$p_net_capture_z           <- tmp$par$p_net_capture_z

save(sra_ini, file = file.path(dat_dir, "sra_ini.rda"))

# Plot fit of initial values ----

# Observed vs expected captures
cdat1 <- data.frame(obs = sra_dat$live_captures_i, var_id = 1:sra_dat$n_i, species = sra_dat$species_i)
cdat2 <- data.frame(obs = sra_dat$dead_captures_i, var_id = 1:sra_dat$n_i, species = sra_dat$species_i)

clive1 <- tmp$par$mu_live_captures_i %>%
  reshape2::melt(varnames = "var_id") %>%
  mutate(grp = "Alive")
cdead1 <- tmp$par$mu_dead_captures_i %>%
  reshape2::melt(varnames = "var_id") %>%
  mutate(grp = "Dead")

clive2 <- clive1 %>%
  left_join(cdat1, by = "var_id") %>%
  group_by(species, grp) %>%
  summarise(obs = sum(obs), 
            hat = sum(value))

cdead2 <- cdead1 %>%
  left_join(cdat2, by = "var_id") %>%
  group_by(species, grp) %>%
  summarise(obs = sum(obs), 
            hat = sum(value))

df1 <- rbind(clive2, cdead2) %>%
  ungroup()

stopifnot(sum(df1$obs) == sum(cdat1$obs) + sum(cdat2$obs))

#gg <- ggplot(data = df1, aes(x = obs, y = hat)) +
#  geom_point(size = 3) +
#  geom_abline(slope = 1) +
#  facet_wrap(grp ~ ., ncol = 1) +
#  theme_bw() +
#  labs(x = "Observed number of captures", y = "Predicted number of captures", colour = NULL) +
#  theme(aspect.ratio = 1)
#ggsave(gg, filename = "sra_ini_fit.png")

rm(cdead2, cdead1, clive2, clive1, cdat1, cdat2)


gg <- ggplot(data = df1, aes(x = obs, y = hat)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  facet_wrap(grp ~ ., ncol = 1) +
  theme_bw() +
  labs(x = "Observed number of captures", y = "Predicted number of captures", colour = NULL) +
  theme(aspect.ratio = 1)
ggsave(gg, file = file.path(dat_dir, "fit_initial.png"))

# Print initial deaths and captures
emp <- data.frame(id = sra_dat$species_i, captures = sra_dat$captures_i, dead_captures = sra_dat$dead_captures_i) %>% group_by(id) %>% summarise(observed_captures = sum(captures), observed_dead_captures = sum(dead_captures))
tab <- data.frame(id = 1:71, species = species_groups$code_sra, predicted_captures = round(tmp$par$observed_captures_s), predicted_dead_captures = round(tmp$par$observed_dead_captures_s), risk = format(signif(tmp$par$risk_ratio_s, 2), nsmall = 2)) %>% left_join(emp, by = "id") %>% dplyr::select(species, observed_captures, predicted_captures, observed_dead_captures, predicted_dead_captures, risk)
write.csv(tab, file = file.path(dat_dir, "sra_ini_risk.csv"), row.names = FALSE)



