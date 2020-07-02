library(ProjectTemplate)
load.project()
theme_set(theme_bw())

library(brms)
library(bayesplot)
library(tidybayes)
options(mc.cores=parallel::detectCores())
bname<-tools::file_path_sans_ext(basename(this.file.name()))
#stop()

# enabls --force 
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options)
  uncache.all(base = bname)


cbind(baseline, part="baseline") %>% 
  rbind(cbind(stimulation,part="stimulation")) -> d

d.nback=get.nback(d, nback=25, which.apen.m = 2)

#========================
## function applied to every model for fitting and plotting
#========================
fit_and_plot <- function(mod.name,frm,load.only=F,plot.only.new=T,init="random"){
  #mod.name = formula.name.fname(frm)
  is.new=TRUE
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = d.nback, family =cumulative("probit"), init=init) %>%
      add_loo() %>% add_waic()
    assign(mod.name, mod, envir=.GlobalEnv)
    cache.var(mod.name, bname)
  } else {
    mod <- load.cache.var(mod.name,bname)
    is.new=FALSE
  }
  if(!load.only & ((is.new & plot.only.new) | (!plot.only.new))  ){
  pdf(plot.filename(sprintf("diag_%s.pdf", mod.name),bname), width=5, height=5)
  mcmc_rhat(brms::rhat(mod)) %>% print
  mcmc_neff(brms::neff_ratio(mod)) %>% print
  dev.off()
  
  
  mcmc_intervals_data(as.matrix(mod), prob_outer = 0.95) %>%
    filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
    ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
    geom_pointrange(position=position_dodge(width=0.2))+
    coord_flip()
  ggsave(plot.filename(sprintf("coef_%s.pdf",mod.name),bname), width=9,height=6)
  
  fit=mod
  nrep=100
  pred=predict(fit)
  
  d.nback %>%
    cbind(
      replicate(n=nrep, apply(pred, 1, function(x){sample(1:4,1, prob=x)})) 
    )   %>%
    gather(sim.n,sim.response, 13:(13+nrep-1)) %>%
    group_by(condition,part, sim.n) %>%
    do({
      tibble(response=1:4,n=tabulate(.$sim.response, nbins=4))
    }) -> d.nback.pred
  
  d.tab=d.nback %>% group_by(condition,part) %>%
    do({
      v=as.numeric(data.frame(.)[,"probe.response"])
      tibble(response=1:4,n=tabulate(v, nbins=4))
    })

  d.nback.pred %>% ungroup %>% 
    ggplot(aes(x=factor(response),y=n,color=condition))+
    geom_bar(data=d.tab, mapping=aes(fill=condition), stat="identity",position = position_dodge(width=1), alpha=0.2)+
    #geom_violin(aes(group=interaction(stim_setting,response),color=NULL),fill="grey",color=0, alpha=1, position=position_dodge(width=1))+
    stat_summary(fun.data = mean_qi,  position=position_dodge(width=1), geom="pointrange")+
    #facet_wrap(~question,ncol=1) +
    facet_grid(part~.) +
    labs(x="Response",y="Number of subjects",
         title=sprintf("%s: Posterior predictive", fit$formula$resp), 
         subtitle=toString(capture.output(fit$formula)))
  
  ggsave(plot.filename(sprintf("ppred_%s.pdf",mod.name),bname), width=9,height=6)
  }
  return(mod)
}


#========================
## model definitions
#========================
models <- list(
  formula(probe.response ~ (1|subj/part)),
  formula(probe.response ~ zbv+(1|subj/part)),
  formula(probe.response ~ zlog.apen+(1|subj/part)),
  formula(probe.response ~ part+(1|subj/part)),
  formula(probe.response ~ zbv+zlog.apen+(1|subj/part)),
  formula(probe.response ~ zbv+zlog.apen+part+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+(1|subj/part)),                                             # 09
  formula(probe.response ~ part+condition+condition:part+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zbv*condition+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zlog.apen*condition+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zbv*zlog.apen*condition+(1|subj/part)),                     # 13
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zbv*condition*part+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zlog.apen*condition*part+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+condition+condition:part+zlog.apen*condition*part+zbv*condition*part+(1|subj/part)), # 16
  formula(probe.response ~ part+probeix+condition+condition:part+(1|subj/part)),                                                   # 17
  formula(probe.response ~ part+probeix+condition+condition:part+condition:probeix+(1|subj/part)),                             
  formula(probe.response ~ part+probeix+condition+condition:part+condition:probeix+condition:part:probeix+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+probeix+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+probeix+condition+(1|subj/part)),
  formula(probe.response ~ zbv*zlog.apen+part+probeix+condition+condition:part+(1|subj/part))
  #formula(probe.response ~ zbv*zlog.apen+part+condition:part+zlog.apen:part+zbv:part+zbv:part:condition + zlog.apen:part:condition+(1|subj/part)),  ## doesn't fit
  #formula(probe.response ~ zbv*zlog.apen+part+condition:part+zlog.apen:part+zbv:part+zbv:part:condition + zlog.apen:part:condition+zlog.apen:zbv:part:condition+(1|subj/part))  ## doesn't fit
)

descriptions=c("Null", "BV", "AE", "part", "BV + AE", "BV + AE + part", "BV x AE", "BV x AE + part",
               "BV x AE + part + stim", "BV x AE + part x stim", "part x stim", 
               "BV x AE + part x stim + BV x stim","BV x AE + part x stim + AE x stim",
               "BV x AE + part x stim + BV x AE x stim", "BV x AE + part x stim + BV x part x stim",
               "BV x AE + part x stim + AE x part x stim", "BV x AE + part x stim + BV x part x stim + AE x part x stim",
               "part x stim + trial", "part x stim + trial x stim", "part x stim x trial", "BV x AE + part + trial",
               "BV x AE + part + stim + trial", "BV x AE + part x stim + trial")
#description=c("Null", "Lab", "Stim(ulation)", "Imp(edance)", "Stim+Lab", "Stim+Imp", "Lab+Imp", "Lab+Stim+Imp", "Stim$\\times$Lab+Imp", "Stim$\\times$Lab$\\times$Imp",  "Stim+Lab", "Stim+Lab")

names(models) <- sprintf("mod%02i", 0:(length(models)-1))

#========================
## fit models
#========================
#library(parallel)
#fit_and_plot("mod0", models[[1]])

library(pbapply)
models.wrap <- map2(names(models), models, ~ list(mod.name=.x, mod=.y))
models.fitted=pblapply(models.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod)}, cl=10)

#mod19=fit_and_plot("mod19", models$mod19, init=0)
#mod22=fit_and_plot("mod22", models$mod22, init=0)

#models.fitted=map2(names(models), models, fit_and_plot)
#models.fitted=future_map2(names(models), models, fit_and_plot, .progress = T)
names(models.fitted) <- names(models)

#map_df(models.fitted, nsamples)
#mod=models.fitted[["mod20"]]
#mod=update(mod, iter=2000, chains=4)
#mod20=mod
#cache.var("mod20", bname)
#models.fitted[["mod20"]]=mod
#uncache.var("mod20", bname)


#models.fitted=map2(names(models), models, fit_and_plot, load.only=T)
#models.fitted=mcmapply(fit_and_plot, names(models), models, mc.cores = 12, SIMPLIFY = FALSE)

#stop()
#========================
## model-selection
#========================

#uncache.var("loos", base=bname)
#uncache.var("mod.weights", base=bname)
loos=if.cached.load("loos",
                    invoke(loo_wrapper, .x = models.fitted, model_names = names(models.fitted)),
                    base=bname)


# Bayesian R2
r2s=pblapply(models.fitted, bayes_R2, cl=22)


#map(map(models.fitted, ~ .x$loo), dim)
#r=loo::loo_model_weights(x=map(models.fitted, ~ .x$loo))
stop()
mod.weights = if.cached.load("mod.weights",
                             map_df(c("loo", "waic", "loo2"), function(strat) {
                               r = invoke(
                                 model_weights_wrapper,
                                 .x = models.fitted,
                                 weights = strat,
                                 model_names = names(models.fitted)
                               )
                               bind_cols(strategy = strat, data.frame(t(r)))
                             }), bname)

sink(report.filename("modsel.log", bname))

print(loos)
print(loo::compare(x=loos$loos))
as.data.frame(loos$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod.weights)
sink(NULL)



## model - weight plot

mod.desc=data.frame(mod=names(models.fitted), descriptions)
map_df(c("loo","loo2","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions) %>%
  #mutate(loo2="",waic="") %>%
  gather(strategy,descriptions,loo,loo2,waic) -> mod.desc

mod.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","loo2")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="loo2" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy!="WAIC") %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

ggsave(plot.filename("study3_model_weights.pdf", bname), width=10,height=5)

mod.weights[,-1] %>% as.matrix %>% t %>% data.frame %>% 
  setNames(mod.weights$strategy) %>%
  rownames_to_column() -> modw.df
modw.df %>% arrange(desc(loo)) %>% head(2)
modw.df %>% arrange(desc(loo2)) %>% head(2)

## LOOIC
loo::compare(x = loos$loos) %>% data.frame %>% rownames_to_column(var = "mod") %>%
  arrange(mod) %>% 
  mutate(rel.looic=looic-first(looic)) %>%
  ggplot(aes(mod,rel.looic))+
  geom_bar(stat="identity", fill="grey")+
  geom_text(
    aes(fill = NULL, label = frm),
    y = 0,
    data = data.frame(mod = names(models), frm =
                        map_chr(
                          models, ~ paste(format(.x), sep = "", collapse = "")
                        )),
    hjust = "right"
  ) +
  #geom_pointrange(aes(ymin=rel.looic-se_looic,ymax=rel.looic+se_looic))+
  coord_flip()

ggsave(plot.filename("rel_looic.pdf", bname), width=12,height=9)



#========================
## inspection of the model
#========================
stop()
test.models = c("mod10", "mod09", "mod13", "mod11", "mod12", "mod14", "mod15", "mod16")
map2_df(test.models,
        models.fitted[test.models],
        ~ bind_cols(
          mod = .x,
          hypothesis(.y, "partstimulation:conditionreal<0", alpha = .05)$hypothesis
        ))


## all models that have the interaction
map2_df(names(models.fitted),
        models.fitted,
        function(modname,mod){
          if("partstimulation:conditionreal" %in% row.names(fixef(mod))){
            return(bind_cols(mod=modname,hypothesis(mod, "partstimulation:conditionreal<0", alpha = .05)$hypothesis))
          } else{
            return(NULL);
          }
        }) -> d.ia
d.ia

d.ia$Evid.Ratio %>% summary

summary(mod22)$fixed %>% data.frame %>% rownames_to_column() %>%
  mutate(summary=sprintf("$b=%.2f\\ [%.2f, %.2f]$", Estimate, l.95..CI, u.95..CI)) %>%
  select(rowname,summary)

hypothesis(mod22e, c("zbv>0", "zlog.apen<0", "partstimulation>0", "probeix>0", 
                     "conditionreal>0", "zbv:zlog.apen>0", "partstimulation:conditionreal<0"))

pred=with(models.fitted, predict(mod22))
library(ggforce)
d.nback %>% bind_cols(data.frame(pred) %>% setNames(c("pprob1","pprob2","pprob3","pprob4"))) %>% 
  #filter(subj==4) %>%
  arrange(part,probeix) %>%
  mutate(probeix=if_else(part=="stimulation", probeix+9, probeix)) %>%
  group_by(condition,part,probeix) %>% 
  mutate(cond.subj=1:n()) %>% ungroup %>%
  arrange(subj) -> d.tmp


d.tmp %>%
  ggplot(aes(probeix, probe.response, color = part)) +
  geom_point(aes(y = ppred, size = probpred, color=NULL),
             color="grey", alpha=0.2,
             data = d.tmp %>%
               gather(ppred, probpred, starts_with("pprob")) %>%
               separate(ppred,into = c("blub","ppred"), sep = 5) %>%
               mutate(ppred=as.integer(ppred))
             )+
  geom_point(aes(y = bestpred),
             color="orange", alpha=0.2, size=10,
             data = bind_cols(d.tmp, bestpred=apply(d.tmp[13:16], 1, which.max))
  )+
  geom_point()+theme_bw() -> p

npages=n_pages(p+facet_grid_paginate(cond.subj~condition, nrow=5, ncol=2,page=1))


pdf(plot.filename("ppred_subj_mod22.pdf", bname))
map(1:npages,
    ~ print(
      p + facet_grid_paginate(
        cond.subj ~ condition,
        nrow = 5,
        ncol = 2,
        page = .x
      )
    ))
dev.off()

####
mod=models.fitted[["mod22"]]
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9 ) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
  geom_pointrange(position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red")+
  labs(y="Coefficient")

## nya SfN
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9 ) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  mutate(parameter=fct_relevel(parameter, "b_zbv:zlog.apen", after=5)) %>%
  mutate(parameter=fct_recode(parameter, 
                              `Variability`="b_zbv",
                              `Entropy`="b_zlog.apen",
                              `Block`="b_partstimulation",
                              `Trial`="b_probeix",
                              `Stimulation`="b_conditionreal",
                              `Variability x Entropy`="b_zbv:zlog.apen",
                              `Block x Stimulation`="b_partstimulation:conditionreal")) %>%
  ggplot(aes(y=m, x=parameter))+
  geom_pointrange(aes(ymin=l,ymax=h), position=position_dodge(width=0.2), color="black", size=2,fatten=0.8)+
  geom_pointrange(aes(ymin=ll,ymax=hh), position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red", size=2, alpha=0.2)+
  labs(y="Coefficient",x="Predictor")+
  #annotate("text", y=-0.5, x=7.25, label="P(b<0)=0.97, Evidence Ratio=30.7",hjust = 0)+
  theme(axis.text.y = element_text(angle=0, hjust=1, size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16)) -> p1
p1
ggsave("graphs/blockxstim.pdf", width=7, height=4)



## coefficients from all relative models
map_df(c("zbv", "zlog.apen", "partstimulation:conditionreal"), function(var) {
  map2_df(names(models.fitted), models.fitted, function(mod.name, mod) {
    var = sprintf("b_%s", var)
    if (var %in% get_variables(mod)) {
      v = as.matrix(mod)[, var]
      data.frame(
        mod = mod.name,
        mean = mean(v),
        lower = hdi(v)[1],
        upper = hdi(v)[2]
      )
    } else {
      data.frame(mod = mod.name)
    }
  }) %>%
    cbind(var = var)
}) %>%
  mutate(best.mod = case_when(mod == "mod22" ~ "LOO",
                              T ~ "none")) %>%
  mutate(var=fct_recode(var, `Block x Stimulation`="partstimulation:conditionreal",
                        `Variability`="zbv",
                        `Entropy`="zlog.apen")) %>%
  mutate(best.mod=ordered(best.mod, c("LOO","none"))) %>%
  mutate(won=if_else(best.mod=="none", F,T)) %>%
  ggplot(aes(x=var,y=mean,color=best.mod)) + 
  geom_pointrange(aes(ymin=lower,ymax=upper,group=mod,size=won), 
                  position=position_dodge(width=0.4), data=(. %>% filter(!won)), alpha=0.5)+
  geom_pointrange(aes(ymin=lower,ymax=upper,group=mod,size=won), alpha=1, 
                  position=position_dodge(width=0.4), data=(. %>% filter(won)))+
  scale_color_manual(values=c("orange", "lightblue", "pink"))+
  scale_size_manual(values=c(0.5,1))+
  coord_flip()+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(y="Coefficient", x="", color="Selected Model")+guides(size=F)+
  theme(legend.position = c(.2,.8),
        legend.background = element_blank())+
  theme(axis.text.y = element_text(angle=0, hjust=1, size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16)) -> p2#-> p.mod.coef
p2
ggsave("graphs/blockxstim_all.pdf", width=7, height=4)

#mod22e = update(mod22, iter=8000, init=0)


library(patchwork)
p=(p1+labs(tag="A")+p2+labs(tag="B"))
ggsave(plot=p, filename = plot.filename("study3_coeffs.pdf",bname), width=12, height=4)


## ntrials between probes

d %>% filter(stimulus=="probe1") %>%
  group_by(subj,part) %>%
  mutate(ntrials=diff(c(1,trial))) %>%
  summarize(min(ntrials),max(ntrials),mean(ntrials))

