
# Dragonfly ??
get_fishing_group_prev <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA, iwl=NA) {
    ifelse(fishery %in% 'FLAT', 20,     ## Flatfish trawl
    ifelse(fishery %in% 'INST' | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S'), 1, ## Inshore trawl
    ifelse(class == 'S' & fishery %in% 'BNSB', 4, ## Small vessel bluenose bottom longline
    ifelse(class == 'S' & fishery %in% c('HAPB', 'MINB'), 5, ## Small vessel bottom longline targetting hapuka and minor species
    ifelse(class == 'S' & fishery %in% c('LINB'), 23,    ## Small vessel bottom longline targetting ling
    ifelse(class == 'S' & fishery %in% c('SNAB'), 6, ## Small vessel snapper
    ifelse(class == 'L' & method == 'BLL', 9, ## Autoline
    ifelse((class == 'L' | (is.na(class) & fishery %in% 'STNS')) & method == 'SLL', 10, ## Large vessel surface longline
    ifelse(class == 'S' & method == 'SLL' & fishery != 'SWOS', 11, ## Small vessel surface longline targetting tuna and other minor species
    ifelse(class == 'S' & method == 'SLL' & fishery %in% 'SWOS', 22, ## Small vessel surface longline targetting swordfish
    ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(15169,12487,12906,12903,15289,20687,6620,11036,15042,15398,13006,6061,15585,11644,8591,20876,5995,6610,15014,15500,6618,6138,15039,13313,5530,5933,6154,21482) | (vessel == 12368 & as.Date(start_date) >= as.Date('2001-01-01'))), 12, ## Meal capable middle-depths trawl
    ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(3725, 15532, 8609,8601,333,6129,11138,5250,1193,1195,8800,327,13106,360,3763,9259,12600,11338,5262,5247,8700,20809,359,804) | (vessel == 6096 & as.Date(start_date) < as.Date('2007-01-01'))), 13, ## Fresher middle depths
    ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(3704, 5558, 8040, 21356, 20610,20884,20864,6984,5663,13521,6944,5458,6473,6645,5981,8804,5906,5921,5803,13009,1356,12599,6730,6489,3351,20766,5605,1282,15256,9921) | (vessel == 12368 & as.Date(start_date) < as.Date('2001-01-01')) | (vessel == 6096 & as.Date(start_date) >= as.Date('2007-01-01'))), 2, ## Freezer middle depths
    ifelse((is.na(class) | class == 'L') & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT'), 3, ## Unclassified middle-depths
    ifelse(fishery %in% 'SBWT', 15, ## Southern blue whiting
    ifelse(fishery %in% 'SCIT', 16, ## Scampi
    ifelse(fishery %in% 'MACT', 17, ## Mackerel
    ifelse(fishery %in% 'SQUT', 18, ## Squid
    ifelse(fishery %in% 'DPWT', 19, ## Deepwater
    ifelse(method == 'SN', 21, ## Set net
    0))))))))))))))))))))
}

# Webber 2020
get_fishing_group_new <- function(method = NA, fishery = NA, target = NA, vessel = NA, class = NA, start_date = NA, length = NA) {
    ### TRAWL
    # ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length < 17, 24, ## Inshore trawl + FLAT < 17 m
    # ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length >= 17, 25, ## Inshore trawl + FLAT >= 17 & < 28 m
    ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('00-06', '06-17'), 24, ## Inshore trawl + FLAT < 17 m
    ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('17-28', '28-43', '43+'), 25, ## Inshore trawl + FLAT >= 17 & < 28 m
    ifelse(fishery %in% 'SBWT', 15, ## Southern blue whiting
    ifelse(fishery %in% 'SCIT', 16, ## Scampi
    ifelse(fishery %in% 'MACT', 17, ## Mackerel
    ifelse(fishery %in% 'SQUT', 18, ## Squid
    ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(1282,1356,3351,3704,5458,5530,5558,5605,5663,5803,5906,5921,5933,5981,5995,6061,6138,6154,6473,6489,6610,6618,6620,6645,6730,6944,6984,8040,8591,8804,9921,11036,11644,12368,12487,12599,12903,12906,13006,13009,13313,13521,15014,15039,15042,15169,15256,15289,15398,15500,15585,20610,20687,20766,20864,20876,20884,21356,21482) | (vessel == 6096 & as.Date(start_date) >= as.Date('2007-01-01'))), 26, ## Non-fresher middle-depths trawl
    ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(3725, 15532, 8609,8601,333,6129,11138,5250,1193,1195,8800,327,13106,360,3763,9259,12600,11338,5262,5247,8700,20809,359,804) | (vessel == 6096 & as.Date(start_date) < as.Date('2007-01-01'))), 13, ## Fresher middle depths
    ifelse(fishery %in% 'DPWT', 19, ## Deepwater
    ifelse((is.na(class) | class == 'L') & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT'), 3, ## Unclassified middle-depths
    #### BLL
    ifelse(class == 'S' & fishery %in% 'BNSB', 4, ## Small vessel bluenose bottom longline
    ifelse(class == 'S' & fishery %in% 'SNAB', 6, ## Small vessel snapper
    ifelse(class == 'S' & method == 'BLL' & target %in% c('LIN', 'RIB'), 27, ## Small vessel bottom longline targetting ling or ribaldo
    ifelse(class == 'S' & fishery %in% c('HAPB', 'MINB'), 28, ## Other small vessel bottom longline
    ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(3761,5541,5947,6976,12595,13315,20914,20648,8601,21514,20809,20914) | (vessel == 15168 & as.Date(start_date) < as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) < as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) < as.Date('2006-01-01'))), 29, ## Large BLL without IWL
    ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(20670,21072) | (vessel == 15168 & as.Date(start_date) >= as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) >= as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) >= as.Date('2006-01-01'))), 30, ## Large BLL with IWL
    ### SLL
    ifelse(class == 'S' & method == 'SLL' & fishery %in% 'SWOS', 22, ## Small vessel surface longline targetting swordfish
    ifelse(class == 'S' & method == 'SLL' & fishery != 'SWOS', 11, ## Small vessel surface longline targetting tuna and other minor species
    ifelse((class == 'L' | (is.na(class) & fishery %in% 'STNS')) & method == 'SLL', 10, ## Large vessel surface longline
    ### SN
    ifelse(method %in% 'SN', 21, ## Set net
    0))))))))))))))))))))
}

# Same as Webber 2020 but assigns id_fishing_method not fishing group code
get_fishing_group_base <- function(method = NA, fishery = NA, target = NA, vessel = NA, class = NA, start_date = NA, length = NA) {

### TRAWL
# ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length < 17, 24, ## Inshore trawl + FLAT < 17 m
# ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length >= 17, 25, ## Inshore trawl + FLAT >= 17 & < 28 m
ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('00-06', '06-17'), 11, #24, ## Inshore trawl + FLAT < 17 m
ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('17-28', '28-43', '43+'), 12, #25, ## Inshore trawl + FLAT >= 17 & < 28 m
ifelse(fishery %in% 'SBWT', 13, #15, ## Southern blue whiting
ifelse(fishery %in% 'SCIT', 14, #16, ## Scampi
ifelse(fishery %in% 'MACT', 15, #17, ## Mackerel
ifelse(fishery %in% 'SQUT', 16, #18, ## Squid
ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(1282,1356,3351,3704,5458,5530,5558,5605,5663,5803,5906,5921,5933,5981,5995,6061,6138,6154,6473,6489,6610,6618,6620,6645,6730,6944,6984,8040,8591,8804,9921,11036,11644,12368,12487,12599,12903,12906,13006,13009,13313,13521,15014,15039,15042,15169,15256,15289,15398,15500,15585,20610,20687,20766,20864,20876,20884,21356,21482) | (vessel == 6096 & as.Date(start_date) >= as.Date('2007-01-01'))), 17, #26, ## Non-fresher middle-depths trawl
ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(3725, 15532, 8609,8601,333,6129,11138,5250,1193,1195,8800,327,13106,360,3763,9259,12600,11338,5262,5247,8700,20809,359,804) | (vessel == 6096 & as.Date(start_date) < as.Date('2007-01-01'))), 18, #13, ## Fresher middle depths
ifelse(fishery %in% 'DPWT', 19, #19, ## Deepwater
ifelse((is.na(class) | class == 'L') & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT'), 20, #3, ## Unclassified middle-depths
#### BLL
ifelse(class == 'S' & fishery %in% 'BNSB', 1, #4, ## Small vessel bluenose bottom longline
ifelse(class == 'S' & fishery %in% 'SNAB', 2, #6, ## Small vessel snapper
ifelse(class == 'S' & method == 'BLL' & target %in% c('LIN', 'RIB'), 3, #27, ## Small vessel bottom longline targetting ling or ribaldo
ifelse(class == 'S' & fishery %in% c('HAPB', 'MINB'), 4, #28, ## Other small vessel bottom longline
ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(3761,5541,5947,6976,12595,13315,20914,20648,8601,21514,20809,20914) | (vessel == 15168 & as.Date(start_date) < as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) < as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) < as.Date('2006-01-01'))), 6, #29, ## Large BLL without IWL
ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(20670,21072) | (vessel == 15168 & as.Date(start_date) >= as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) >= as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) >= as.Date('2006-01-01'))), 5, #30, ## Large BLL with IWL
### SLL
ifelse(class == 'S' & method == 'SLL' & fishery %in% 'SWOS', 7, #22, ## Small vessel surface longline targetting swordfish
ifelse(class == 'S' & method == 'SLL' & fishery != 'SWOS', 8, #11, ## Small vessel surface longline targetting tuna and other minor species
ifelse((class == 'L' | (is.na(class) & fishery %in% 'STNS')) & method == 'SLL', 9, #10, ## Large vessel surface longline
### SN
ifelse(method %in% 'SN', 10, #21, ## Set net
0))))))))))))))))))))
}


bll_groups_original <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA) {
    ifelse(class == 'S' & fishery %in% 'BNSB', 4, ## Small vessel bluenose bottom longline
           ifelse(class == 'S' & fishery %in% 'SNAB', 6, ## Small vessel snapper
                  ifelse(class == 'S' & method == 'BLL' & target %in% c('LIN', 'RIB'), 27, ## Small vessel bottom longline targetting ling or ribaldo
                         ifelse(class == 'S' & fishery %in% c('HAPB', 'MINB'), 28, ## Other small vessel bottom longline
                                ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(3761,5541,5947,6976,12595,13315,20914,20648,8601,21514,20809,20914) | (vessel == 15168 & as.Date(start_date) < as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) < as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) < as.Date('2006-01-01'))), 29, ## Large BLL without IWL
                                       ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(20670,21072) | (vessel == 15168 & as.Date(start_date) >= as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) >= as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) >= as.Date('2006-01-01'))), 30, ## Large BLL with IWL
                                              ifelse(method == 'BLL', 0, # unclassified BLL
                                                     NA)))))))
}

## Autoliners
bll_autoline_test <- function(method, vessel, date) {
    vessel_autoline <- c(1699, 3761, 5541, 5892, 5947, 6952, 6976, 12595, 13315, 15165, 15168,
                         15229, 15544, 20662, 20670, 20809, 20914, 21014, 21182, 21514)
    method %in% 'BLL' & (
        vessel %in% vessel_autoline |
            ## Vessels that moved from manual to autoline
            (vessel %in% 15229 & as.Date(date) >= as.Date('2014-08-18')) |
            (vessel %in% 21182 & as.Date(date) >= as.Date('2012-08-24')) |
            ## Vessels that moved from autoline to manual
            (vessel %in% 20546 & as.Date(date) < as.Date('2013-07-04')) |
            (vessel %in% 8601 & as.Date(date) < as.Date('2019-01-01')))
}

## Manual BLL
bll_manual_test <- function(method, vessel, date) {
    vessel_manual <- c(26, 34, 57, 71, 83, 163, 168, 174, 176, 177, 226, 251, 316, 319, 324,
                       337, 346, 409, 439, 452, 474, 529, 539, 548, 580, 584, 586, 760, 770,
                       788, 794, 805, 807, 842, 848, 856, 861, 1132, 1140, 1279, 1382, 1426,
                       1484, 1493, 1541, 1558, 1602, 1642, 1760, 1772, 1829, 1842, 1861, 1865,
                       1971, 1972, 2015, 2089, 2125, 2127, 2143, 2159, 2226, 2618, 2634, 2770,
                       2775, 2795, 2987, 2992, 2998, 3043, 3205, 3472, 3491, 3554, 3586, 3604,
                       3748, 3777, 3811, 3823, 3892, 4077, 4115, 4574, 4594, 4672, 4854, 5114,
                       5347, 5463, 5606, 6040, 6128, 6179, 6216, 6998, 8612, 9038, 15131,
                       15258, 15317, 15511, 15572, 15647, 20769, 20983, 21059, 21348, 21427,
                       21480, 21599)
    
    method %in% 'BLL' & (
        vessel %in% vessel_manual |
            ## Vessels that moved from manual to autoline
            (vessel %in% 15229 & as.Date(date) < as.Date('2014-08-18')) |
            (vessel %in% 21182 & as.Date(date) < as.Date('2012-08-24')) |
            ## Vessels that moved from autoline to manual
            (vessel %in% 20546 & as.Date(date) >= as.Date('2013-07-04')) |
            (vessel %in% 8601 & as.Date(date) >= as.Date('2019-01-01')))
}


bll_groups_option_II <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA) {
    ## Small vessel autoliner
    ## LIN & RIB
    ifelse(class == 'S' & method %in% 'BLL' & bll_autoline_test(method, vessel, start_date) &
               target %in% c('LIN', 'RIB'), 911,
           
           ## others
           ifelse(class == 'S' & method %in% 'BLL' & bll_autoline_test(method, vessel, start_date) &
                      !target %in% c('LIN', 'RIB'), 912,
                  
                  ## Small vessel manual
                  ## light gear
                  ifelse(class == 'S' & method %in% 'BLL' & bll_manual_test(method, vessel, start_date) &
                             target %in% c('SNA', 'BCO', 'GUR', 'KAH', 'RRC', 'RSN', 'SPO', 'SUN', 'TAR',
                                           'TRE', 'KIN', 'ELE', 'JDO', 'SNX'), 901,
                         
                         ## heavy gear - LIN & RIB
                         ifelse(class == 'S' & method %in% 'BLL' & bll_manual_test(method, vessel, start_date) &
                                    target %in% c('LIN', 'RIB'), 903,
                                
                                ## heavy gear - others
                                ifelse(class == 'S' & method %in% 'BLL' & bll_manual_test(method, vessel, start_date) &
                                           target %in% c('BAS', 'BNS', 'HAP', 'HPB', 'SCH',
                                                         'BYX', 'SPE', 'ALB', 'OSD', 'SND', 'SPD'), 902,
                                       
                                       ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(3761,5541,5947,6976,12595,13315,20914,20648,8601,21514,20809,20914) | (vessel == 15168 & as.Date(start_date) < as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) < as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) < as.Date('2006-01-01'))), 929, ## Large BLL without IWL
                                              ifelse(class == 'L' & method == 'BLL' & (vessel %in% c(20670,21072) | (vessel == 15168 & as.Date(start_date) >= as.Date('2003-03-01')) | (vessel == 20662 & as.Date(start_date) >= as.Date('2003-11-01')) | (vessel == 15165 & as.Date(start_date) >= as.Date('2006-01-01'))), 930, ## Large BLL with IWL
                                                     ifelse(method == 'BLL', 900, # unclassified BLL
                                                            NA))))))))
}

## Names of fisheries groups
lk_bll_groups_option_II <- data.frame(
    fishing_group = c(929, 930, 911, 912, 901, 902, 903, 900),
    fishing_group_name = c('BLL-large-auto-without-IWL', 'BLL-large-auto-with-IWL',
                           'BLL-small-auto-ling-ribaldo', 'BLL-small-auto-other',
                           'BLL-small-manual-light', 'BLL-small-manual-heavy',
                           'BLL-small-manual-ling-ribaldo',
                           'BLL-unclassified'))


trawl_groups_original <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA) {
    ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('00-06', '06-17'), 24, ## Inshore trawl + FLAT < 17 m
           ifelse((fishery %in% c('INST','FLAT') | (fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & class == 'S')) & length %in% c('17-28', '28-43', '43+'), 25, ## Inshore trawl + FLAT >= 17 & < 28 m
                  ifelse(fishery %in% 'SBWT', 15, ## Southern blue whiting
                         ifelse(fishery %in% 'SCIT', 16, ## Scampi
                                ifelse(fishery %in% 'MACT', 17, ## Mackerel
                                       ifelse(fishery %in% 'SQUT', 18, ## Squid
                                              ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(1282,1356,3351,3704,5458,5530,5558,5605,5663,5803,5906,5921,5933,5981,5995,6061,6138,6154,6473,6489,6610,6618,6620,6645,6730,6944,6984,8040,8591,8804,9921,11036,11644,12368,12487,12599,12903,12906,13006,13009,13313,13521,15014,15039,15042,15169,15256,15289,15398,15500,15585,20610,20687,20766,20864,20876,20884,21356,21482) | (vessel == 6096 & as.Date(start_date) >= as.Date('2007-01-01'))), 26, ## Non-fresher middle-depths trawl
                                                     ifelse(class == 'L' & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT') & (vessel %in% c(3725, 15532, 8609,8601,333,6129,11138,5250,1193,1195,8800,327,13106,360,3763,9259,12600,11338,5262,5247,8700,20809,359,804) | (vessel == 6096 & as.Date(start_date) < as.Date('2007-01-01'))), 13, ## Fresher middle depths
                                                            ifelse(fishery %in% 'DPWT', 19, ## Deepwater
                                                                   ifelse((is.na(class) | class == 'L') & fishery %in% c('LINT', 'MIDT', 'HAKT', 'HOKT'), 3, ## Unclassified middle-depths
                                                                          ifelse(method == 'Trawl', 0, # other unclassified trawls
                                                                                 NA)))))))))))
}

## Dave's option 2A - squid and mackerel fishery groups not defined spatially
trawl_groups_option_IIA <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA, freezer=NA) {
    deepwater_spp <- c('BOE', 'CDL', 'OEO', 'ORH', 'SSO')
    
    ifelse(length %in% c('00-06', '06-17') & !target %in% c('SCI', deepwater_spp), 824, ## Trawl < 17 m (not scampi or deep-water)
           ifelse(length %in% '17-28' & !target %in% c('SCI', deepwater_spp), 825, ## Trawl 17 - 28 m (not scampi or deep-water)
                  ifelse(length %in% c('28-43', '43+') & !freezer & !target %in% c('SCI', deepwater_spp), 826, ## Trawl 28 m (no freezer, not scampi or deep-water)
                         
                         ifelse(target %in% 'SCI', 816, ## Scampi
                                ifelse(target %in% deepwater_spp, 819, ## Deepwater species
                                       ifelse(target %in% 'SBW' & length %in% c('28-43', '43+'), 815, ## Southern blue whiting
                                              
                                              ifelse(length %in% c('28-43', '43+') & target %in% 'JMA', 817, ## Mackerel
                                                     ifelse(length %in% c('28-43', '43+') & target %in% 'SQU',  818, ## Squid
                                                            
                                                            ifelse(length %in% c('28-43', '43+') & freezer, 830, ## non-fresher everything else
                                                                   
                                                                   ifelse(method == 'Trawl', 800, # other unclassified trawls
                                                                          NA))))))))))
}

## Names of fisheries groups
lk_trawl_groups_option_IIA <- data.frame(
    fishing_group = c(824, 825, 826, 816, 819, 815, 817, 818, 830, 800),
    fishing_group_name = c('Trawl-small-inshore<17m', 'Trawl-small-inshore_17-28m',
                           'Trawl-large-fresher', 'Trawl-scampi', 'Trawl-deepwater',
                           'Trawl-southern-blue-whiting', 'Trawl-mackerel',
                           'Trawl-squid', 'Trawl-large-freezer', 'Trawl-unclassified'))

## Dave's option 2B - squid and mackerel fishery groups defined spatially
trawl_groups_option_IIB <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA, freezer=NA, trawl_area=NA) {
    deepwater_spp <- c('BOE', 'CDL', 'OEO', 'ORH', 'SSO')
    
    ifelse(length %in% c('00-06', '06-17') & !target %in% c('SCI', deepwater_spp), 824, ## Trawl < 17 m (not scampi or deep-water)
           ifelse(length %in% '17-28' & !target %in% c('SCI', deepwater_spp), 825, ## Trawl 17 - 28 m (not scampi or deep-water)
                  ifelse(length %in% c('28-43', '43+') & !freezer & !target %in% c('SCI', deepwater_spp), 826, ## Trawl 28 m (no freezer, not scampi or deep-water)
                         
                         ifelse(target %in% 'SCI', 816, ## Scampi
                                ifelse(target %in% deepwater_spp, 819, ## Deepwater species
                                       ifelse(target %in% 'SBW' & length %in% c('28-43', '43+'), 815, ## Southern blue whiting
                                              
                                              ifelse(length %in% c('28-43', '43+') & target %in% 'JMA' & trawl_area %in% 'WCNI9', 817, ## Mackerel
                                                     ifelse(length %in% c('28-43', '43+') & target %in% 'SQU' & trawl_area %in% c('SQUAK6', 'STEW5'), 818, ## Squid
                                                            
                                                            ifelse(length %in% c('28-43', '43+') & freezer, 830, ## non-fresher everything else
                                                                   
                                                                   ifelse(method == 'Trawl', 800, # other unclassified trawls
                                                                          NA))))))))))
}

## Names of fisheries groups (same as IIA)
lk_trawl_groups_option_IIB <- data.frame(
    fishing_group = c(824, 825, 826, 816, 819, 815, 817, 818, 830, 800),
    fishing_group_name = c('Trawl-small-inshore<17m', 'Trawl-small-inshore_17-28m',
                           'Trawl-large-fresher', 'Trawl-scampi', 'Trawl-deepwater',
                           'Trawl-southern-blue-whiting', 'Trawl-mackerel',
                           'Trawl-squid', 'Trawl-large-freezer', 'Trawl-unclassified'))


## Dave's option 3 - with squid and mackerel fishery groups defined spatially
trawl_groups_option_III <- function(method=NA, fishery=NA, target=NA, vessel=NA, class=NA, start_date=NA, length=NA, freezer=NA, waste_strategy=NA) {
    deepwater_spp <- c('BOE', 'CDL', 'OEO', 'ORH', 'SSO')
    
    ifelse(length %in% c('00-06', '06-17'), 824, ## Trawl < 17 m
           ifelse(length %in% '17-28' & !target %in% c('SCI', deepwater_spp), 825, ## Trawl 17 - 28 m (not scampi or deep-water)
                  ifelse(length %in% c('28-43', '43+') & !freezer & !target %in% c('SCI', deepwater_spp), 826, ## Trawl 28 m (no freezer, not scampi or deep-water)
                         
                         ifelse(target %in% 'SCI', 816, ## Scampi
                                ifelse(!length %in% c('00-06', '06-17') & target %in% deepwater_spp, 819, ## Deepwater species
                                       
                                       ifelse(length %in% c('28-43', '43+') & freezer & waste_strategy %in% 'Meal', 8301, ## Large freezers - meal
                                              ifelse(length %in% c('28-43', '43+') & freezer & waste_strategy %in% 'Batch', 8302, ## Large freezers - batch
                                                     ifelse(length %in% c('28-43', '43+') & freezer & waste_strategy %in% 'Mince', 8303, ## Large freezers - mince
                                                            
                                                            ifelse(method == 'Trawl', 800, # other unclassified trawls
                                                                   NA)))))))))
}

## Names of fisheries groups
lk_trawl_groups_option_III <- data.frame(
    fishing_group = c(824, 825, 826, 816, 819, 8301, 8302, 8303, 800),
    fishing_group_name = c('Trawl-small-inshore<17m', 'Trawl-17m-to-28m',
                           'Trawl-large-fresher', 'Trawl-scampi', 'Trawl-deepwater',
                           'Trawl-large-freezer-meal', 'Trawl-large-freezer-batch',
                           'Trawl-large-freezer-mince', 'Trawl-unclassified'))

# set net
setnet_groups_option_II <- function(method = NA) {
    ifelse(method == 'SN', 700, NA)
}

lk_setnet_groups_option_II <- data.frame(fishing_group = 700, fishing_group_name = 'SN-unclassified')

# PS
purseseine_groups_option_II <- function(method = NA) {
    ifelse(method == 'PS', 600, NA)
}

lk_purseseine_groups_option_II <- data.frame(fishing_group = 600, fishing_group_name = 'PS-unclassified')

# SLL
sll_groups_option_II <- function(class = NA, method = NA) {
    ifelse(class == 'S' & method == 'SLL', 500, ifelse(class == 'L' & method == 'SLL', 550, NA))
}

lk_sll_groups_option_II <- data.frame(fishing_group = c(500, 550), fishing_group_name = c('SLL-small', 'SLL-large'))


# combine names
fishing_groups_option_II <- rbind(lk_trawl_groups_option_IIA, lk_bll_groups_option_II, lk_setnet_groups_option_II, lk_purseseine_groups_option_II, lk_sll_groups_option_II)

fishing_groups_option_II <- dplyr::mutate(fishing_groups_option_II, method = dplyr::case_when(grepl("BLL",   fishing_group_name) ~ "BLL", 
                                                               grepl("SLL",   fishing_group_name) ~ "SLL",
                                                               grepl("PS",    fishing_group_name) ~ "PS",
                                                               grepl("SN",    fishing_group_name) ~ "SN",
                                                               grepl("Trawl", fishing_group_name) ~ "Trawl",
                                                               TRUE ~ "NA"))

fishing_groups_option_II <- dplyr::mutate(fishing_groups_option_II, fishing_group_label = dplyr::case_when(fishing_group == 930 ~ "Large Autoline with IWL",
                                                                                 fishing_group == 929 ~ "Large Autoline",
                                                                                 fishing_group == 911 ~ "Small Autoline (LIN, RIB)",
                                                                                 fishing_group == 912 ~ "Small Autoline",
                                                                                 fishing_group == 902 ~ "Small Manual (heavy)",
                                                                                 fishing_group == 901 ~ "Small Manual (light)",
                                                                                 fishing_group == 903 ~ "Small Manual (LIN, RIB)",
                                                                                 fishing_group == 900 ~ "BLL (unclassified)",
                                                                                 fishing_group == 600 ~ "PS (unclassified)",
                                                                                 fishing_group == 500 ~ "Small SLL (tuna and swordfish)",
                                                                                 fishing_group == 550 ~ "Large SLL",
                                                                                 fishing_group == 700 ~ "SN (unclassified)",
                                                                                 fishing_group == 819 ~ "Deepwater",
                                                                                 fishing_group == 830 ~ "Large Freezer",
                                                                                 fishing_group == 826 ~ "Large Fresher",
                                                                                 fishing_group == 817 ~ "Mackerel",
                                                                                 fishing_group == 816 ~ "Scampi",
                                                                                 fishing_group == 825 ~ "Small inshore (17-28m)",
                                                                                 fishing_group == 824 ~ "Small inshore (<17m)",
                                                                                 fishing_group == 815 ~ "Southern Blue Whiting",
                                                                                 fishing_group == 818 ~ "Squid",
                                                                                 fishing_group == 800 ~ "Trawl (unclassified)"))

fishing_groups_option_II$id_fishing_group_option_II <- fishing_groups_option_II$fishing_group
fishing_groups_option_II$fishing_group <- NULL

fishing_groups_option_II <- fishing_groups_option_II[order(fishing_groups_option_II$fishing_group_name),]


