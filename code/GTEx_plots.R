library(here)
library(tidyverse)

dfs <- readRDS(here('output', 'gtex_dfs.rds'))



dfs[[1]] %>%
  mutate(setting = 'Setting 1: Normal-Normal') %>%
  ggplot(aes(x=time, y=log(diff_from_maxelbo), color=method))+
  theme_bw()+
  geom_line()+geom_point()+
  facet_wrap(~setting, scales='free')+
  xlab('Time (s)')+ylab('Log difference from maximum ELBO')+
  ggtitle('GTEx: Log Difference from maximum ELBO')
ggsave(here('figures', 'gtex_n-n.pdf'), width=5, height=4)



rbind(mutate(dfs[[2]], setting='Setting 2: PointNormal-PointNormal'),
      mutate(dfs[[3]], setting='Setting 3: PointNormal-PointExponential')) %>%
  ggplot(aes(x=time, y=log(diff_from_maxelbo), col=method))+
  theme_bw()+
  geom_line()+
  facet_wrap(~setting, scales='free')+
  #scale_x_log10()+
  xlab('Time (s)') + ylab('Log difference from maximum ELBO')+
  ggtitle('GTEx: Log Difference from maximum ELBO')
ggsave(here('figures', 'gtex_pn-pn_pn-pe.pdf'), width=8, height=6)

