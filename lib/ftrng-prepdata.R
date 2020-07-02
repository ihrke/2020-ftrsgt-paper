get.nback.pilot2 <- function(d, nback=20, which.apen.m=3, on.task.crit=2){

  d %>% filter(stimulus=="probe1") %>% 
    mutate(attention=if_else(as.integer(response)<=on.task.crit, "on-task", "off-task")) %>%
    group_by(subj) %>% 
    do({
      df=.
      df %>% mutate(probeix=1:n()) -> df
      dd=d %>% filter(subj==df$subj[1])
      probeix=which(dd$stimulus=="probe1")
      nprobes=length(probeix)
      dd$probeix=rep(c(1:nprobes, -1), diff(c(0,probeix,dim(dd)[1])))
      
      trials=map(1:nback, function(x) df$trial-x) %>% unlist
      dd %>%
        filter(trial %in% trials) %>%
        mutate(focus=map_chr(probeix, function(t){
            df$attention[t]
          }),
          probe.response=map_int(probeix, function(t){
            as.integer(df$response[t])
          })+1)
    }) %>% ungroup -> d.nback
  
  
  d.nback %<>% group_by(subj) %>% 
    filter(stimulus=="tap") %>% mutate(tap=case_when(response=="rctrl" ~ 1,
                                                     response=="lalt" ~ 0,
                                                     TRUE ~ -1)) %>% 
    filter(tap>=0) %>% mutate(tap=as.integer(tap)) %>% ungroup 
  
  
  d.nback = map_df(which.apen.m, ~bind_cols(d.nback, m=rep(.x, dim(d.nback)[1])))

  d.nback %>% 
    group_by(subj, focus, probeix,probe.response,m) %>%
    summarize(apen=apen_int(tap,3)[first(m)+1]) %>%
    mutate(log.apen=log(log(2)-apen)) %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,focus,m,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% 
    group_by(m) %>%
    mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task")) %>%
    ungroup -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}



get.nback <- function(d, nback=20, which.apen.m=3, on.task.crit=1){
  if( !("part" %in% names(d))){
    d$part=1 ## this is for pilot2
  }
  d %>% filter(stimulus=="probe1") %>% 
    mutate(attention=if_else(as.integer(response)<=on.task.crit, "on-task", "off-task")) %>%
    group_by(subj,part) %>% 
    do({
      df=.
      df %>% mutate(probeix=1:n()) -> df
      dd=d %>% filter(subj==df$subj[1], part==df$part[1])
      probeix=which(dd$stimulus=="probe1")
      nprobes=length(probeix)
      dd$probeix=rep(c(1:nprobes, -1), diff(c(0,probeix,dim(dd)[1])))
      
      trials=map(1:nback, function(x) df$trial-x) %>% unlist
      dd %>%
        filter(trial %in% trials) %>%
        mutate(focus=map_chr(probeix, function(t){
          df$attention[t]}),
          probe.response=map_int(probeix, function(t){
            as.integer(df$response[t])
          })+1)
    }) %>% ungroup -> d.nback
  
  
  d.nback %<>% group_by(part, subj) %>% 
    filter(stimulus=="tap") %>% mutate(tap=case_when(response=="rctrl" ~ 1,
                                                     response=="lctrl" ~ 0,
                                                     TRUE ~ -1)) %>% 
    filter(tap>=0) %>% mutate(tap=as.integer(tap)) %>% ungroup 
  
  
  d.nback %>% 
    group_by(part, subj, focus, probeix,probe.response) %>%
    summarize(apen=apen_int(tap,3)[which.apen.m+1]) %>%
    mutate(log.apen=log(log(2)-apen)) %>% left_join(groups,by="subj") %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,part,focus,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>%
    left_join(groups,by=c("subj")) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task"),
           condition=fct_relevel(condition,"sham")) -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}


#
# trial-wise re-arrangement of data.
# use as: 
#
# pilot1 %>% mutate(ISI=as.factor(ISI)) %>% group_by(subj, ISI) %>% 
#   do( rearrange.df(.) ) %>% ungroup %>% 
#   mutate(reltime=resp_time-stim_time) %>% data.frame -> pilot1.bytrial
#
#
rearrange.ftrngt.bytrial.df.onesubj <- function(df){
  stim_times=with(df, time[stimulus==1])
  resp_times=with(df, time[stimulus==0])
  resp_trial=unlist(lapply(resp_times, function(x){ which.min(abs(stim_times+.1-x))}))
  d<-data.frame(
    subj=df$subj[1],
    ISI=df$ISI[1],
    trial=1:length(stim_times),
    stim_time=stim_times,
    resp_time=NA,
    response=NA,
    nresponses=0
  )
  d<-within(d,{
    resp_time[resp_trial]=resp_times;
    response[resp_trial]=with(df, as.character(response[stimulus==0]));
    nresponses=unlist(lapply(1:dim(d)[1], function(i){sum(i==resp_trial)}));
  })
}
