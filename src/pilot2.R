library(ProjectTemplate)
load.project()
theme_set(theme_bw())
library(rstan)
library(tidybayes)
library(ggridges)
library(patchwork)
library(brms)
library(magrittr)
library(broom) 
library(modelr) 
library(purrr)
library(pROC)
library(caret)
library(lmerTest)

options(mc.cores=parallel::detectCores())

bname<-tools::file_path_sans_ext(basename(this.file.name()))


maxM = 5
stop()

# look for missing trials/unusual stuff
pilot2 %>% group_by(subj,stimulus) %>%
  summarise(n=n()) %>%
  spread(stimulus, n) %>% data.frame

## distance bw probes
pilot2 %>% filter(stimulus=="probe1") %>%
  group_by(subj) %>%
  mutate(probediff=diff(c(0,trial))) %>%
  summarise(min(probediff),mean(probediff),max(probediff)) %>%
  data.frame

## all responses
pilot2 %>% filter(stimulus == "tap") %>%
  mutate(response.int =
           as.integer(as.factor(as.character(response))) - 1) %>%
  mutate(subj = as.factor(subj)) %>%
  group_by(subj) %>%
  do(data.frame(apen = apen_int(.$response.int, maxM), m = as.factor(0:maxM))) %>%
  filter(m != 0) %>%
  mutate(logapen=-log(log(2)-apen)) %>% ungroup %>% droplevels-> d

d %>% gather(which.apen, apen, apen, logapen) %>%
  ggplot(aes(x=apen,y=m))+
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE)+
  facet_wrap(~which.apen, scales="free")


d %>% select(-apen) %>% mutate(m=sprintf("m%i",m)) %>% spread(m,logapen) %>%
  select(-subj) %>%
  GGally::ggpairs()

d.nback=get.nback.pilot2(pilot2, nback=20, which.apen=2)
d.nback %>% gather(var,val,zlog.apen,zbv) %>%
  ggplot(aes(x=focus,y=val, color=var))+
  stat_summary(fun.data=mean_se, geom="pointrange")+
  stat_summary(fun.y=mean, mapping=aes(group=var), fun.args=list(na.rm=T), geom="line")


library(BayesFactor)
d.nback %>% select(subj,probeix,focus,zlog.apen, zbv) %>% group_by(subj,focus) %>%
  summarize(apen=mean(zlog.apen, na.rm=T), bv=mean(zbv, na.rm=T)) %>% gather(var,val,apen,bv) %>%
  unite(var_focus, var,focus) %>% spread(var_focus, val) %>% as.data.frame %>% na.omit -> d.bf
write_csv(d.bf, path="data/export/pilot2_bf.csv")
anovaBF(val ~ var*focus, whichRandom="subj", data=d.bf)
d.nback$m

## all m
######################
d=map_df(1:3, ~ cbind(m=.x,get.nback.pilot2(pilot2, nback=20, which.apen=.x)))

d %>% mutate(zbv=ifelse(m==1, zbv, NA)) %>%
  gather(var, val, zlog.apen, zbv) %>%
  mutate(m=as.factor(m)) %>%
  mutate(var=fct_recode(var,BV="zbv",AE="zlog.apen")) %>%
  ggplot(aes(x = focus, y = val, color = var)) +
  stat_summary(
    fun.data = mean_se,
    mapping = aes(group = interaction(var, m)),
    position = position_dodge(width = 0.3),
    geom = "pointrange"
  ) +
  stat_summary(
    fun.y = mean,
    mapping = aes(group = interaction(var, m), linetype = m),
    fun.args = list(na.rm = T),
    position = position_dodge(width = 0.3),
    geom = "line"
  ) +
  labs(x="Focus",y="Z-score",color="Variable",linetype="AE (m)")+
  coord_cartesian(expand=F, xlim = c(0.8,2.5))+
  theme(legend.position = c(0.9,0.5),
        legend.background = element_blank())-> p.cross.m
p.cross.m

## all nback
######################
d=map_df(seq(10,25,by=5), ~ cbind(nback=.x,get.nback.pilot2(pilot2, nback=.x, which.apen=2)))

d %>% gather(var, val, zlog.apen, zbv) %>%
  mutate(nback=factor(nback)) %>%
  mutate(var=fct_recode(var,BV="zbv",`AE (m=2)`="zlog.apen")) %>%
  ggplot(aes(x = focus, y = val, color = var)) +
  stat_summary(
    fun.data = mean_se,
    mapping = aes(group = interaction(var, nback)),
    position = position_dodge(width = 0.1),
    geom = "pointrange"
  ) +
  stat_summary(
    fun.y = mean,
    mapping = aes(group = interaction(var, nback), linetype = nback),
    fun.args = list(na.rm = T),
    position = position_dodge(width = 0.1),
    geom = "line"
  ) +
  coord_cartesian(expand=F, xlim = c(0.8,2.5))+
  guides(color=guide_legend(order=1))+
  labs(x="Focus",y="Z-score",linetype="nback",color="Variable")+
  theme(legend.position = c(0.88,0.5),
        legend.background = element_blank())->  p.cross.nback
p.cross.nback
#========================
## function applied to every model for fitting and plotting
#========================
fit_and_plot <- function(mod.name,frm,data=NULL, load.only=F,plot.only.new=T){
  #mod.name = formula.name.fname(frm)
  is.new=TRUE
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = data, family =cumulative("probit")) %>%
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
  }
  return(mod)
}


#========================
## model definitions
#========================
models <- list(
  formula(probe.response ~ (1|subj)),
  formula(probe.response ~ zprobeix + (1|subj)),
  formula(probe.response ~ zbv+(1|subj)),
  formula(probe.response ~ zlogapen1+(1|subj)),
  formula(probe.response ~ zlogapen2+(1|subj)),
  formula(probe.response ~ zlogapen1+zlogapen2+(1|subj)),
  formula(probe.response ~ zbv+zlogapen1+(1|subj)),
  formula(probe.response ~ zbv+zlogapen2+(1|subj)),
  formula(probe.response ~ zbv+zlogapen1+zlogapen2+(1|subj)),
  formula(probe.response ~ zbv+zlogapen1+zprobeix+(1|subj)),
  formula(probe.response ~ zbv+zlogapen2+zprobeix+(1|subj)),
  formula(probe.response ~ zbv+zlogapen1+zlogapen2+zprobeix+(1|subj)),
  formula(probe.response ~ zbv*zlogapen1+(1|subj)),
  formula(probe.response ~ zbv*zlogapen2+(1|subj)),
  formula(probe.response ~ zbv*zlogapen1*zlogapen2+(1|subj)),
  formula(probe.response ~ zbv*zlogapen1+zprobeix+(1|subj)),
  formula(probe.response ~ zbv*zlogapen2+zprobeix+(1|subj)),
  formula(probe.response ~ zbv*zlogapen1*zlogapen2+zprobeix+(1|subj))
)
descriptions = c(
  "Null",
  "Trial",
  "BV",
  "AE (m=1)",
  "AE (m=2)",
  "AE(m=1,2)",
  "BV+AE(m=1)",
  "BV+AE(m=2)",
  "BV+AE(m=1,2)",
  "BV+AE(m=1)+Trial",
  "BV+AE(m=2)+Trial",
  "BV+AE(m=1,2)+Trial",
  "BV x AE(m=1)",
  "BV x AE(m=2)",
  "BVx AE(m=1,2)",
  "BV x AE(m=1)+Trial",
  "BV x AE(m=2)+Trial"
  #"BV x AE(m=1,2) x Trial"
)

names(models) <- sprintf("mod%02i", 0:(length(models)-1))

#========================
## fit models
#========================
d.nback=get.nback.pilot2(pilot2, nback=25, which.apen.m=c(1,2))
d.nback = within(d.nback, {
  zlog.apen[is.na(zlog.apen)]=0;
})
d.nback %<>% 
  mutate(zprobeix=(probeix-mean(probeix))/sd(probeix)) %>%
  select(-log.apen, -apen) %>%
  mutate(m=sprintf("zlogapen%i",m)) %>%
  spread(m,zlog.apen) 


library(furrr)
plan(multiprocess, workers = 12)
models.fitted=map2(names(models), models, data=d.nback, fit_and_plot, load.only=T)
#models.fitted=future_map2(names(models), models, fit_and_plot, data=d.nback, .progress = T)
names(models.fitted) <- names(models)
models.fitted=models.fitted[-length(models.fitted)]

## check description
map2_df(models.fitted,
        descriptions,
        ~ data.frame(desc = .y, form = reduce(deparse(.x$formula$formula), paste)))

#========================
## model-selection
#========================

#uncache.var("loos", base=bname)
#uncache.var("mod.weights", base=bname)
loos=if.cached.load("loos",
                    invoke(loo_wrapper, .x = models.fitted, model_names = names(models.fitted)),
                    base=bname)

#r=loo::loo_model_weights(x=map(models.fitted, ~ .x$loo))

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
print(loo::compare(x=loos[-length(loos)]))
as.data.frame(loos$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod.weights)
sink(NULL)

mod.desc=data.frame(mod=names(models.fitted), descriptions)
map_df(c("loo","loo2","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions) %>%
  #mutate(loo2="",waic="") %>%
  gather(strategy,descriptions,loo,loo2,waic) -> mod.desc

## model - weight plot
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
        strip.placement = "inside") -> p.mod.prob
p.mod.prob
ggsave(plot.filename("model_weights.pdf", bname), width=12,height=9)
 
mod.weights[,-1] %>% as.matrix %>% t %>% data.frame() %>% set_colnames(mod.weights$strategy)%>% rownames_to_column()-> modw.df
modw.df %>% arrange(desc(loo)) %>% head(2)
modw.df %>% arrange(desc(loo2)) %>% head(2)

## winning model
#-----------------------
mod10=models.fitted[["mod10"]]
mod16=models.fitted[["mod16"]]
modm10=as.array(mod10)
modm16=as.array(mod16)
library(bayesplot)
mcmc_areas_ridges(modm10, regex_pars = "b_.*")
mcmc_areas_ridges(modm16, regex_pars = "b_.*")

hypothesis(mod10, c("zbv>0", "zlogapen2<0"))
hypothesis(mod16, c("zbv>0", "zlogapen2<0", "zprobeix>0", "zbv:zlogapen2>0"))

## coefficients from all relative models
map_df(c("zbv", "zlogapen1", "zlogapen2"), function(var) {
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
  mutate(best.mod = case_when(mod == "mod10" ~ "LOO",
                              mod == "mod16" ~ "pseudo-BMA",
                              T ~ "none")) %>%
  mutate(best.mod=ordered(best.mod, c("LOO","pseudo-BMA","none")),
         var=fct_recode(var, `AE(m=1)`="zlogapen1",
                        `AE(m=2)`="zlogapen2",
                        `BV`="zbv")) %>%
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
  theme(legend.position = c(.8,.8),
        legend.background = element_blank()) -> p.mod.coef




library(gridExtra)
library(grid)
arrangeGrob(p.cross.m+labs(tag="A"), 
             p.cross.nback+labs(tag="B"), 
             p.mod.prob+labs(tag="C"), 
             p.mod.coef+labs(tag="D"), 
             layout_matrix=matrix(c(1,3,2,4),ncol=2))->p
ggsave(filename = plot.filename("pilot2_results.pdf",bname), plot = p, width=22, height=18, units = "cm")

arrangeGrob(p.cross.m+labs(tag="A"), 
            p.cross.nback+labs(tag="B"), 
            layout_matrix=matrix(c(1,2),ncol=2))->p
ggsave(filename = plot.filename("pilot2_cross.pdf",bname), plot = p, width=22, height=9, units = "cm")



arrangeGrob(p.mod.prob+labs(tag="A"), 
            p.mod.coef+labs(tag="B"), 
            layout_matrix=matrix(c(1,2),ncol=2))->p
ggsave(filename = plot.filename("pilot2_model_weights.pdf",bname), plot = p, width=22, height=9, units = "cm")





## LOOIC
loo::compare(x = loos[-length(loos)]) %>% data.frame %>% rownames_to_column(var = "mod") %>%
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
