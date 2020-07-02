library(ProjectTemplate)
load.project()
theme_set(theme_bw())
library(rstan)
library(tidybayes)
library(ggridges)
library(patchwork)
library(brms)
library(magrittr)
options(mc.cores=parallel::detectCores())

bname<-tools::file_path_sans_ext(basename(this.file.name()))

maxM = 3
stop()

# num trials per ISI
pilot1 %>% filter(stimulus==1) %>%
  group_by(subj,ISI) %>%
  summarise(n=n()) %>% group_by(ISI) %>%
  summarise(mean(n))

## -> problem: much more trials in ISI=0.3 condition 
## -> solution 1: get only first 250 trials
## -> solution 2: bootstrap randomly pick 250 trials from each ISI's seq and average apen

## all responses
pilot1 %>% filter(stimulus == 0) %>%
  mutate(response.int =
           as.integer(as.factor(as.character(response))) - 1) %>%
  mutate(subj = as.factor(subj)) %>%
  group_by(subj, ISI) %>%
  do(data.frame(apen = apen_int(.$response.int, maxM), m = as.factor(0:maxM))) %>%
  filter(m != 0) %>%
  mutate(logapen=-log(log(2)-apen)) -> d


### -----------------------------------
## distribution of the ApEns
### -----------------------------------


d %>% gather(which.apen, apen, apen, logapen) %>%
  droplevels() %>%
  group_by(which.apen, m) %>%
  nest() %>%
  mutate(shapiro=map(data, ~ broom::tidy(shapiro.test(.$apen)))) %>%
  unnest(shapiro) %>%
  select(which.apen, m, W=statistic, p=p.value) %>%
  mutate(summary=sprintf("W=%.2f, p=%s", W,pformat(p)),
         x=if_else(which.apen=="apen", 0.0, 0))-> d.apen.shapiro


d %>% gather(which.apen, apen, apen, logapen) %>%
  ggplot(aes(apen,y=m,fill=m))+
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE)+
  geom_text(mapping=aes(x=x, y=as.integer(as.factor(m))-0.1, label=summary), data=d.apen.shapiro, size=4, hjust="left")+
  facet_wrap(~which.apen,scales="free",ncol=1)+
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )-> p.apen.distr


d %>% gather(which.apen, apen, apen, logapen) %>%
  ggplot(aes(sample=apen, color=m))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~which.apen,scales="free",ncol=1)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
  )-> p.apen.qq

p=p.apen.distr+labs(title="Distribution of ApEn and -log(log(2)-ApEn)")+p.apen.qq+labs(caption=bname)
ggsave(plot.filename("apen_distr.pdf",bname), plot=p+labs(caption=""), width=7,height=4)


### -----------------------------------
## performance by ISI
### -----------------------------------


d %>% 
  gather(var,value,apen,logapen) %>%
  ggplot(aes(ISI, value)) +
  geom_hline(yintercept = log(2), color = 'grey') +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  geom_line(aes(color = subj)) + 
  facet_grid(var ~ m, scales="free")+
  labs(title="Raw and log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: different amount of trials in each ISI condition
decrease for longer ISIs probably because less trials (more trials, easier to get high apen)")
ggsave(plot.filename("apen_by_isi_ind.pdf",bname), width=9,height=6)

d %>% 
  gather(var,value,apen,logapen) %>%
  ggplot(aes(ISI, value)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") + 
  stat_summary(fun.y = mean, geom = "line") + 
  facet_grid(var ~ m, scales="free")+
  labs(title="Raw and log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: different amount of trials in each ISI condition
decrease for longer ISIs probably because less trials (more trials, easier to get high apen)")
ggsave(plot.filename("apen_by_isi.pdf",bname), width=9,height=6)


d %>% 
  ggplot(aes(ISI, logapen,color=m)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") + 
  stat_summary(fun.y = mean, geom = "line") + 
  labs(title="Log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: different amount of trials in each ISI condition
decrease for longer ISIs probably because less trials (more trials, easier to get high apen)")
ggsave(plot.filename("logapen_by_isi_color.pdf",bname), width=6,height=4)



### -----------------------------------
## ApEn/ISI compared to self-eval
### -----------------------------------
pilot1 %>% filter(stimulus==2) %>% mutate(subj=as.factor(subj)) %>%
  mutate(evaluation=(as.integer(as.character(response))+1)) %>%
  ggplot(aes(ISI,evaluation))+
  #geom_line(aes(color=subj))+
  stat_summary(fun.data=mean_se, geom="pointrange")+
  stat_summary(fun.y=mean, geom="line")-> p.eval.isi

ggsave(plot.filename("evaluation_by_isi.pdf",bname), width=6,height=4)

#ylim(0,1)
pilot1 %>% filter(stimulus==2) %>% mutate(subj=as.factor(subj)) %>%
  mutate(evaluation=(as.integer(as.character(response))+1)) %>%
  select(subj,ISI,evaluation) %>%
  full_join(d) %>%
  ggplot(aes(y=evaluation,x=logapen,color=m))+
  geom_point(position=position_jitter(height=0.1),alpha=0.3)+
  geom_smooth(method="lm", se=F,fullrange=T)+
  facet_grid(.~ISI)+
  coord_cartesian(ylim=c(1,5))

ggsave(plot.filename("cor_logapen_evaluation_by_isi.pdf",bname), width=9,height=4)


pilot1 %>% filter(stimulus==2) %>% mutate(subj=as.factor(subj)) %>%
  mutate(evaluation=(as.integer(as.character(response))+1)) -> d.tmp

library(brms)
mod=brm(evaluation ~ factor(ISI)+(1|subj), data=d.tmp)
hypothesis(mod, c("factorISI0.75>factorISI0.5", 
                  "factorISI1>factorISI0.75", 
                  "factorISI1.25>factorISI1",
                  "factorISI0.75>factorISI0.75"
                  ))
mod2=brm(evaluation ~ as.numeric(ISI)+(1|subj), data=d.tmp)
hypothesis(mod2, "as.numericISI>0")

### -----------------------------------
## -> solution 1: get only first 250 trials
### -----------------------------------

pilot1 %>% filter(stimulus == 0) %>% 
  mutate(response.int =
           as.integer(as.factor(as.character(response))) - 1) %>%
  mutate(subj = as.factor(subj)) %>%
  group_by(subj, ISI) %>%
  mutate(trial=1:n()) %>% 
  ungroup %>%
  filter(trial<=250) %>%
  group_by(subj,ISI) %>%
  do(data.frame(apen = apen_int(.$response.int, maxM), m = as.factor(0:maxM))) %>%
  filter(m != 0) %>%
  mutate(logapen=-log(log(2)-apen)) -> d

d %>% 
  gather(var,value,apen,logapen) %>%
  ggplot(aes(ISI, value)) +
  geom_hline(yintercept = log(2), color = 'grey') +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  geom_line(aes(color = subj)) + 
  facet_grid(var ~ m, scales="free")+
  labs(title="Raw and log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: only the first 250 trials have been used for each ISI condition")
ggsave(plot.filename("apen_by_isi_ind_first250.pdf",bname), width=9,height=6)

d %>% 
  gather(var,value,apen,logapen) %>%
  ggplot(aes(ISI, value)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") + 
  stat_summary(fun.y = mean, geom = "line") + 
  facet_grid(var ~ m, scales="free")+
  labs(title="Raw and log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: Note: only the first 250 trials have been used for each ISI condition")
ggsave(plot.filename("apen_by_isi_first250.pdf",bname), width=9,height=6)

d %>% 
  ggplot(aes(ISI, logapen,color=m)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") + 
  stat_summary(fun.y = mean, geom = "line") + 
  labs(title="Log-transformed ApEn as a function of ISI", 
       subtitle=bname,
       caption="Note: different amount of trials in each ISI condition
decrease for longer ISIs probably because less trials (more trials, easier to get high apen)")-> p.apen.isi.250

ggsave(plot=p.apen.isi.250, filename = plot.filename("logapen_by_isi_color_first250.pdf",bname), width=6,height=4)

#ylim(0,1)
pilot1 %>% filter(stimulus==2) %>% mutate(subj=as.factor(subj)) %>%
  mutate(evaluation=(as.integer(as.character(response))+1)) %>%
  select(subj,ISI,evaluation) %>%
  full_join(d) -> d.eval

d.eval %>%
  ggplot(aes(y=evaluation,x=logapen,color=m))+
  geom_point(position=position_jitter(height=0.1),alpha=0.3)+
  geom_smooth(method="lm", se=F,fullrange=T)+
  facet_grid(.~ISI)+
  coord_cartesian(ylim=c(1,5))

ggsave(plot.filename("cor_logapen_evaluation_by_isi.pdf",bname), width=9,height=4)

### -----------------------------------
# bayesian correlations per ISI
mod.cor=stan_model(file="lib/robust_correlation.stan")
d.eval %>% 
  select(m,ISI,logapen,evaluation) %>%
  group_by(m,ISI) %>%
  do({
    mat=as.matrix(.[,3:4])
    fit=sampling(mod.cor, data=list(x=mat,N=nrow(mat)), iter=2000, warmup=1000, chains=4)    
    rho=as.vector(rstan::extract(fit, pars="rho")[[1]] )
    data.frame(
      mean_rho=mean(rho),
      lower=hdi(rho)[1],
      upper=hdi(rho)[2]
    )
  }) -> d.eval.cor
  
d.eval.cor %>%
  ggplot(aes(ISI,mean_rho,ymin=lower,ymax=upper,color=m))+
  geom_pointrange(position=position_dodge(width=0.1))+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(y="Correlation",
       title="Correlation between measured and self-perceived Randomness",
       caption=bname) -> p.corr.eval.apen


ggsave(p.corr.eval.apen, plot.filename("cor_random_self.pdf",bname), width=6,height=4)

### -----------------------------------
## SD(ITI) as function of ISI
### -----------------------------------

pilot1 %>% group_by(subj,ISI) %>%
  do({
    dx=data.frame(.)
    dx %>% filter(stimulus==0) %>% mutate(dtime=c(0,diff(time))) %>% filter(dtime>0) %>% pull(dtime) -> x
    
    data.frame(
      mean=mean(x),sd=sd(x),cv=sd(x)/mean(x),min=min(x), max=max(x)
    )
    }) %>% ungroup -> d.sd

## distribution

d.sd %>%
  nest(-ISI) %>%
  mutate(shapiro=map(data, ~ broom::tidy(shapiro.test(.$cv)))) %>%
  unnest(shapiro) %>%
  select(ISI, W=statistic, p=p.value) %>%
  mutate(summary=sprintf("W=%.2f, p=%s", W,pformat(p)))-> d.sd.shapiro

d.sd %>% mutate(ISI=factor(ISI)) %>%
  ggplot(aes(x=cv,y=ISI))+
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE)+
  geom_text(mapping=aes(x=0.25, y=as.integer(as.factor(ISI))+0.5, label=summary), data=d.sd.shapiro, size=4, hjust="left")+
  labs(y="ISI",x="CV(ITI)") -> p.sd.dist

d.sd %>% mutate(ISI=factor(ISI)) %>%
  ggplot(aes(sample=cv,color=ISI))+
  geom_qq()+
  geom_qq_line()+
  theme(legend.position = c(0.2,0.6),
        legend.background = element_blank())-> p.sd.qq
p=p.sd.dist+p.sd.qq
ggsave(plot.filename("cv_iti_distr.pdf",bname), plot=p, width=6,height=3)

## performance
d.sd %>%
  ggplot(aes(x=ISI, y=cv))+
  geom_line(aes(group=subj),color="lightblue")+
  stat_summary(fun.data=mean_se, geom="pointrange")+
  stat_summary(fun.y=mean, geom="line")+
  coord_cartesian(ylim=c(0.05,0.3))+
  labs(y="CV (ITI)") -> p.sd
p.sd

ggsave(plot.filename("cv_isi.pdf",bname), plot=p.sd, width=4,height=3)

(p.sd.dist+p.sd.qq+ p.sd + plot_layout(ncol=3) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))) -> p.sd.paper

sc=1.3
ggsave(plot=p.sd.paper, filename=plot.filename("cv_paper.pdf", bname), units = "cm", width=sc*21, height=sc*6)



### -----------------------------------
## models logapen - ISI  for detecting optimal ISI
### -----------------------------------

mod00=if.cached.load("mod00", 
                     brm(logapen ~  (1|subj), data=d), 
                     base=bname)

mod0a=if.cached.load("mod0a", 
                    brm(logapen ~ poly(ISI,1)  + (1|subj), data=d), 
                    base=bname)
mod0b=if.cached.load("mod0b", 
                     brm(logapen ~ m  + (1|subj), data=d), 
                     base=bname)
mod0c=if.cached.load("mod0c", 
                     brm(logapen ~ poly(ISI,1) + m + (1|subj), data=d), 
                     base=bname)

mod1=if.cached.load("mod1", 
                    brm(logapen ~ poly(ISI,1) * m + (1|subj), data=d), 
                    base=bname)
mod2=if.cached.load("mod2", 
                    brm(logapen ~ poly(ISI,2) * m + (1|subj), data=d),
                    base=bname)
mod2b=if.cached.load("mod2b", 
                     brm(logapen ~ ISI * m + I(ISI^2)*m + (1|subj), data=d), ## just a raw 2nd order poly model (poly() uses orthogonalized polys)
                     base=bname)
mod3=if.cached.load("mod3", brm(logapen ~ poly(ISI,3) * m + (1|subj), data=d),
                    base=bname)

mod00 %<>% add_loo()
mod0a %<>% add_loo()
mod0b %<>% add_loo()
mod0c %<>% add_loo()
mod1 %<>% add_loo()
mod2 %<>% add_loo()
mod2b %<>% add_loo()
mod3 %<>% add_loo()
descriptions=c(mod00="Null", mod0a="ISI", mod0b="m", mod0c="ISI + m", mod1="ISI x m", mod2="quadratic(ISI) x m", mod3="cubic(ISI) x m")


loos=loo(mod00, mod0a, mod0b, mod0c, mod1,mod2b,mod3)
# LOOIC    SE
# mod1        832.71 38.40
# mod2        817.87 37.34
# mod3        823.15 37.27
# mod1 - mod2  14.84 11.33
# mod1 - mod3   9.56 11.83
# mod2 - mod3  -5.28  2.58
weights=model_weights(mod1,mod2b,mod3)
#         mod1         mod2         mod3 
# 2.406100e-01 7.593897e-01 2.110041e-07 

## LOOIC
library(latex2exp)
loo::compare(x = loos$loos) %>% data.frame %>% rownames_to_column(var = "mod") %>%
  arrange(mod) %>% 
  mutate(rel.looic=looic-first(looic)) %>%
  mutate(win=if_else(looic==min(looic), T,F)) %>%
  ggplot(aes(mod,rel.looic))+
  geom_bar(stat="identity", mapping=aes(fill=win))+#"lightblue")+
  geom_text(aes(fill=NULL, label=descriptions), y=0, hjust="right")+
  scale_fill_manual(values=c("lightblue", "orange"))+
  coord_flip()+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank())+labs(x="",y=TeX("$\\Delta$ LOOIC"))-> p.looic


## marginal plots
marginal_effects(mod2b,ask=F)

fe=fixef(mod2b)[,"Estimate"]

# analytical max-point of quadratic model
max.m1=-(fe["ISI"]/(2*fe["IISIE2"]))
max.m2=-( (fe["ISI"]+fe["ISI:m2"])/(2*(fe["IISIE2"]+fe["m2:IISIE2"]) ))
max.m3=-( (fe["ISI"]+fe["ISI:m3"])/(2*(fe["IISIE2"]+fe["m3:IISIE2"]) ))

# analytical maximum of quadratic model for each of the posterior samples

get_variables(mod2b)
mod2b %>% spread_draws(b_ISI, b_IISIE2, `b_ISI:m2`, `b_ISI:m3`, `b_m2:IISIE2`, `b_m3:IISIE2`) %>%
  mutate(max.m1=-b_ISI/(2*b_IISIE2),
         max.m2=-( (b_ISI+`b_ISI:m2`)/(2*(b_IISIE2+`b_m2:IISIE2`))),
         max.m3=-( (b_ISI+`b_ISI:m3`)/(2*(b_IISIE2+`b_m3:IISIE2`)))
         ) %>% 
  gather(m,ISI,starts_with("max")) %>%
  separate(m,c("blub","m"),sep=".m") %>% select(-blub) -> d.mod2b.maxm

d.mod2b.maxm %>%
  ggplot(aes(ISI,y=m))+
  geom_density_ridges()+
  lims(x=c(0.3,1.25))+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  labs(title="Estimated maximum of the ISI-ApEN relationship",
       caption=bname) -> p.mod.max.apenisi
ggsave(plot=p.mod.max.apenisi, filename = plot.filename("mod2b_isi_maxapen.pdf",bname), plot=p1, width=6,height=4)


## which of the measures ISIs is closest?
ISIs=sort(unique(d$ISI))

d.mod2b.maxm %>%
  group_by(m) %>%
  mutate(cISI=factor(ISIs[map_dbl(ISI, ~ closest(ISIs, .x))], levels=ISIs)) %>%
  group_by(m,cISI) %>%
  summarise(n=n()/nsamples(mod2b)) %>%
  ggplot(aes(cISI,n,fill=m))+geom_bar(stat="identity",position = position_dodge())+
  labs(y="Relative Frequency",x="measured ISI",
       title="Closest match from the set of measured ISIs relative to the optimum ISI",
       caption=bname) -> p.mod.closestisi

ggsave(plot=p.mod.closestisi, filename = plot.filename("mod2b_bestisi.pdf",bname),plot=p2, width=6,height=4)


p=p1+p2
ggsave(plot.filename("mod2b_bestisi_joint.pdf",bname),plot=p, width=10,height=4)

### -----------------------------------
## correlations with Baddely 10-digit RNGT
### -----------------------------------
pilot1.rngt %>% group_by(subj) %>% filter(stimulus==0) %>%
  mutate(response.int=as.integer(as.factor(as.character(response)))-1) %>%
  do( data.frame(ApEnRNGT=apen_int(.$response.int, maxM), m=as.factor(0:maxM))) %>% ungroup %>%
  mutate(subj=as.factor(subj)) %>% filter(m!=0) %>%
  mutate(logapen.rngt=-log(log(10)-ApEnRNGT)) %>% full_join(d) -> d.rngt
d.rngt %>%
  group_by(m,ISI) %>%
  do({
    mat=as.matrix(.[,c("logapen","logapen.rngt")])
    fit=sampling(mod.cor, data=list(x=mat,N=nrow(mat)), iter=2000, warmup=1000, chains=4)    
    rho=as.vector(rstan::extract(fit, pars="rho")[[1]] )
    data.frame(
      mean_rho=mean(rho),
      lower=hdi(rho)[1],
      upper=hdi(rho)[2]
    )
  }) -> d.cor.baddeley

#library(gtools)
#map_df(ISIs, ~ cbind(combinations(n = 3, r = 2, repeats.allowed = TRUE), ISI=.x) %>% data.frame) %>%
expand.grid(m1=1:3, m2=1:3, ISIs) %>%
  t %>% as.data.frame %>% future_map_dfr(function(par){
    m.ft=par[1]
    m.bad=par[2]
    cISI=par[3]
    bad <- d.rngt %>% filter(m==m.bad, ISI==cISI) %>% pull(logapen.rngt)
    ft  <- d.rngt %>% filter(m==m.ft,  ISI==cISI) %>% pull(logapen)
    mat=cbind(bad,ft)
    fit=sampling(mod.cor, data=list(x=mat,N=nrow(mat)), iter=2000, warmup=1000, chains=4)    
    rho=as.vector(rstan::extract(fit, pars="rho")[[1]] )
    data.frame(
      m.ft,
      m.bad,
      ISI=cISI,
      mean_rho=mean(rho),
      lower=hdi(rho)[1],
      upper=hdi(rho)[2]
    )
  }) -> d.ccor.baddeley

d.ccor.baddeley %>%
  ggplot(aes(x=ISI,y=mean_rho,ymin=lower,ymax=upper,color=factor(m.ft)))+
  geom_pointrange(position=position_dodge(width=0.15))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_wrap(~m.bad)


d.rngt %>%
  group_by(m,ISI) %>%
  do({
    mat=as.matrix(.[,c("logapen","logapen.rngt")])
    fit=sampling(mod.cor, data=list(x=mat,N=nrow(mat)), iter=2000, warmup=1000, chains=4)    
    rho=as.vector(rstan::extract(fit, pars="rho")[[1]] )
    data.frame(
      mean_rho=mean(rho),
      lower=hdi(rho)[1],
      upper=hdi(rho)[2]
    )
  }) -> d.ccor.baddeley



d.cor.baddeley %>%
  ggplot(aes(ISI,mean_rho,ymin=lower,ymax=upper,color=m))+
  geom_pointrange(position=position_dodge(width=0.1))+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(y="Correlation",
       title="Robust Bayesian correlation between FT-RNGT and Baddeley-task",
       caption=bname) -> p.cor.apen.baddeley
p.cor.apen.baddeley
ggsave(plot=p.cor.apen.baddeley, filename = plot.filename("ft_baddeley_robust.pdf",bname), width=5,height=3)

d.rngt %>%
  group_by(m,ISI) %>%
  do(
    broom::tidy(cor.test(.$logapen,.$logapen.rngt))
  ) %>%
  ggplot(aes(ISI,estimate,color=m,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(position=position_dodge(width=0.1))+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(y="Correlation",
       title="Pearson-correlation between FT-RNGT and Baddeley-task",
       caption=bname)
ggsave(plot.filename("ft_baddeley_pearson.pdf",bname), width=5,height=3)

d.rngt %>%
  ggplot(aes(logapen,logapen.rngt,color=m))+
  geom_point()+
  geom_smooth(method="lm",fullrange=T)+
  facet_grid(m~ISI,scales="free")
ggsave(plot.filename("ft_baddeley_scatter.pdf",bname), width=9,height=5)


d.ccor.baddeley %>%
  mutate(summary=sprintf("ISI=%.02f, $m$=%s: $\\rho=%.2f\\ [%.2f, %.2f]$", ISI, m, mean_rho, lower, upper))


save(p.cor.apen.baddeley, p.looic, p.apen.isi.250, p.eval.isi, p.corr.eval.apen, p.sd, p.apen.distr, p.apen.qq, file = plot.filename("plotobjs.RData", bname))

library(patchwork)
p.apen.distr+labs(tag="A") + theme(plot.tag.position = c(0.05,1))+ 
  (p.apen.qq+labs(tag="B") + theme(plot.tag.position = c(0.05,1)) + theme(legend.position = "none")) + 
  ((p.eval.isi+labs(tag="C"))/
     (p.cor.apen.baddeley+theme(legend.position = "none")+labs(tag="D"))) + 
  ((p.apen.isi.250+labs(title="",subtitle="",caption="")+theme(legend.position = "none")+labs(tag="E"))/
     (p.mod.max.apenisi+labs(title="",subtitle="",caption="")+labs(tag="F"))) + 
  ((p.looic+labs(tag="G"))/
     (p.mod.closestisi+labs(title="",subtitle="",caption="")+labs(tag="H"))) + 
  plot_layout(ncol=5) &     
  theme(plot.margin = unit(c(0,0,0.2,0.4), "cm"),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.tag.position = "topleft", #c(0.5,0.5),
        plot.tag=element_text(face = "bold")) -> p

sc=1.3
ggsave(plot=p, filename=plot.filename("apen.paper.pdf", bname), units = "cm", width=sc*21, height=sc*8)

### TODO: remove distribution info and add cross-correlations with baddeley task
((p.eval.isi+labs(tag="C"))/
    (p.cor.apen.baddeley+theme(legend.position = "none")+labs(tag="D"))) + 
  ((p.apen.isi.250+labs(title="",subtitle="",caption="")+theme(legend.position = "none")+labs(tag="E"))/
     (p.mod.max.apenisi+labs(title="",subtitle="",caption="")+labs(tag="F"))) + 
  ((p.looic+labs(tag="G"))/
     (p.mod.closestisi+labs(title="",subtitle="",caption="")+labs(tag="H")))+
  plot_layout(ncol=3)



library(gridExtra)
library(grid)

theme_multi <- function(...){
  theme(...,
        plot.margin = unit(c(0,0.2,0.1,0.2), "cm"),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.tag.position = "topleft", #c(0.5,0.5),
        plot.tag=element_text(face = "bold"))
}

layoutmat=matrix(c(1,3,5,7,7,
                   8,8,8,7,7,
                   2,4,6,7,7),ncol=5, byrow = T)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


mylegend = g_legend(
  p.apen.isi.250 + labs(colour="AE(m)") + theme_multi(
    legend.direction = "horizontal",
    legend.background = element_rect(
      linetype = "solid",
      color =
        "transparent",
      fill = "transparent"
    )
  )
)

arrangeGrob(p.eval.isi+labs(tag="A")+theme_multi(), 
            p.cor.apen.baddeley+theme_multi(legend.position = "none")+labs(tag="B",y="correlation"), 
            p.apen.isi.250+theme_multi(legend.position = "none", legend.direction="horizontal")+labs(tag="C",y="AE"), 
            p.mod.max.apenisi+labs(tag="E")+theme_multi(), 
            p.looic+labs(tag="D")+theme_multi(), 
            p.mod.closestisi+labs(tag="F")+theme_multi(legend.position="none", 
                                                       legend.background=element_blank(),
                                                       panel.grid=element_blank()), 
            p.sd+labs(tag="G")+theme_multi(),
            mylegend,
            layout_matrix=layoutmat, heights=c(5,1,5))->p
ggsave(filename = plot.filename("pilot1_results.pdf",bname), plot = p, width=25, height=10, units = "cm")

  


###
maxM=3
x=map_df(c(20,50,100,200,250), ~ data.frame(N=.x,t(replicate(n=1000, apen_int(as.integer(runif(.x)>0.5), maxM))))) %>% set_colnames(c("N",seq(0,maxM)))
x %>% gather(m,apen,-N) %>%
  mutate(N=factor(N)) %>%
  mutate(apen=-log(log(2)-apen)) %>%
  ggplot(aes(x=m, y=apen, color=N, group=N))+
  stat_summary(fun.data=mean_qi, geom="pointrange", position=position_dodge(0.3))+
  coord_flip()

x %>%
  mutate(apen=-log(log(2)-apen)) %>%
  ggplot(aes(x=m, y=apen, color=N, group=N))