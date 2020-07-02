library(ProjectTemplate)
load.project()

full_join(demographics,groups) %>%
  mutate(stimulation=as.integer(condition=="real")) %>% 
  select(subj,group.guess=`Antatt gruppe`,stimulation) -> d

chisq.test(d$group.guess, d$stimulation)

library(summarytools)
ctable(d$group.guess, d$stimulation)

library(BayesFactor)
tab=table(d$group.guess, d$stimulation)
1/contingencyTableBF(tab, sampleType = "indepMulti", fixedMargin = "cols", priorConcentration = 1)
