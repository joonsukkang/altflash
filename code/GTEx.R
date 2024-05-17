# run altflash and flashier on GTEx data

library(flashier)
library(here)
library(tidyverse)
source(here('code', 'altflash.R'))
source(here('code', 'solG.R'))
source(here('code', 'solQ.R'))
source(here('code', 'solTAU.R'))
source(here('code', 'altflash-util.R'))


gtex <- readRDS(url('https://github.com/stephenslab/gtexresults/raw/master/data/MatrixEQTLSumStats.Portable.Z.rds'))
gtex <- gtex$strong.z
print(dim(gtex))


fit_df <- function(fit, methodname){
  data.frame(time = fit$vec.time,
             elbo = fit$vec.elbo,
             method = methodname)
}



# Setting 1: n-n
ebnm_fn = ebnm_normal

init <- flash_init(gtex) %>%
  flash_greedy(Kmax = 50L,
               ebnm_fn = ebnm_fn,
               verbose = 0L)
print(init$n_factors)

fit_flashier1 <- flash_backfit_modified(init, extrapolate = TRUE)
fit_flashier2 <- flash_backfit_modified(init, extrapolate = FALSE)

fit_altflash <- flash_init_to_altflash(init,
                                       L_ebnm_fn = ebnm_fn,
                                       F_ebnm_fn = ebnm_fn)
fit_altflash <- altflash_normal(fit_altflash, tol=fit_flashier1$flash_fit$conv.tol/10)

rbind(fit_df(fit_flashier1, 'flashier:extrapolate'),
      fit_df(fit_flashier2, 'flashier:sequential'),
      fit_df(fit_altflash, 'altflash:closedform')) %>%
  mutate(diff_from_maxelbo = max(elbo) - elbo ) -> df1


# Setting 2: pn-pn
ebnm_fn = ebnm_point_normal

init <- flash_init(gtex) %>%
  flash_greedy(Kmax = 50L,
               ebnm_fn = ebnm_fn,
               verbose = 0L)
print(init$n_factors)

fit_flashier1 <- flash_backfit_modified(init, extrapolate = TRUE)
fit_flashier2 <- flash_backfit_modified(init, extrapolate = FALSE)

fit_altflash <- flash_init_to_altflash(init,
                                       L_ebnm_fn = ebnm_fn,
                                       F_ebnm_fn = ebnm_fn)
for (i in 1:5){
  assign(paste0('fit_altflash', i),
         altflash(fit_altflash,
                  solQb_nupdates = i))
}
rbind(fit_df(fit_flashier1, 'flashier:extrapolate'),
      fit_df(fit_flashier2, 'flashier:sequential'),
      fit_df(fit_altflash1, 'altflash:nupdates=1'),
      fit_df(fit_altflash2, 'altflash:nupdates=2'),
      fit_df(fit_altflash3, 'altflash:nupdates=3'),
      fit_df(fit_altflash4, 'altflash:nupdates=4'),
      fit_df(fit_altflash5, 'altflash:nupdates=5')) %>%
  mutate(diff_from_maxelbo = max(elbo) - elbo ) -> df2


# Setting 3: pn-pe
ebnm_fn = list(ebnm_point_normal, ebnm_point_exponential)

init <- flash_init(gtex) %>%
  flash_greedy(Kmax = 50L,
               ebnm_fn = ebnm_fn,
               verbose = 0L)
print(init$n_factors)

fit_flashier1 <- flash_backfit_modified(init, extrapolate = TRUE)
fit_flashier2 <- flash_backfit_modified(init, extrapolate = FALSE)

fit_altflash <- flash_init_to_altflash(init,
                                       L_ebnm_fn = ebnm_fn[[1]],
                                       F_ebnm_fn = ebnm_fn[[2]])
for (i in 1:5){
  assign(paste0('fit_altflash', i),
         altflash(fit_altflash,
                  solQb_nupdates = i))
}
rbind(fit_df(fit_flashier1, 'flashier:extrapolate'),
      fit_df(fit_flashier2, 'flashier:sequential'),
      fit_df(fit_altflash1, 'altflash:nupdates=1'),
      fit_df(fit_altflash2, 'altflash:nupdates=2'),
      fit_df(fit_altflash3, 'altflash:nupdates=3'),
      fit_df(fit_altflash4, 'altflash:nupdates=4'),
      fit_df(fit_altflash5, 'altflash:nupdates=5')) %>%
  mutate(diff_from_maxelbo = max(elbo) - elbo ) -> df3



# save results
saveRDS(list(df1=df1, df2=df2, df3=df3), here('output', 'gtex_dfs.rds'))
