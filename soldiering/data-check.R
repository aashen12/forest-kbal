d1 <- read_csv("blattman.csv")
d2 <- read_csv("blattman_claude.csv")

all(d1$educ == d2$educ)
sum(!is.na(d1$educ))
sum(!is.na(d2$educ))

v1 = d1$log.wage
v2 = d2$lwage_mo



samp <- sample(1:length(v1), size = 25, replace = FALSE)

v1[samp]
v2[samp]

v2[is.na(v2)]
v1[is.na(v2)]

AC_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$", 
                names(d2), value = TRUE)
ACG_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$|^G[1-8]_", 
                 names(d2), value = TRUE)

ctrls <- c("fthr_ed0", "fthr_ed4", "fthr_ed7", "mthr_ed0", "mthr_ed4", "mthr_ed7",
           "no_fthr96", "no_mthr96", "orphan96", "hh_fthr_frm", "hh_size96", 
           "hh_size96_2", "hh_size96_3", "hh_size96_4", "hh_land", "hh_land_2", 
           "hh_land_3", "hh_land_4", "hh_cattle", "hh_cattle_2", "hh_cattle_3", 
           "hh_cattle_4", "hh_stock", "hh_stock_2", "hh_stock_3", "hh_stock_4", "hh_plow")

ctrl_hh <- c("age", "hh_fthr_frm", "hh_size96", "hh_land", "landrich", 
             "hh_cattle", "hh_stock", "hh_plow")

ctrl_i <- c("fthr_ed0", "fthr_ed", "mthr_ed0", "mthr_ed", "no_fthr96", 
            "no_mthr96", "orphan96")

full_covs <- c(AC_vars, ACG_vars, ctrls, ctrl_hh, ctrl_i)

num_vals <- apply(d2 %>% dplyr::select(all_of(full_covs)), 2, function(x) length(unique(x)))
discrete_covs <- names(num_vals[num_vals <= 10])

