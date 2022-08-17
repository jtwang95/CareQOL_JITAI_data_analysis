################################################################
# Codes to performing complete-case analysis in AS project
# Jitao Wang
# 10/30/2021
################################################################

library(tidyverse)
library(ggplot2)
library(geepack)
library(sjPlot)
library(gridExtra)
library(ggpubr)
options(digits=3) 
options(pillar.sigfig=3)

if (.Platform$OS.type == "unix"){
  prefix = "~/Dropbox (University of Michigan)/CareQOL Project/"
} else {
  prefix = "C:/Users/jtwan/Dropbox (University of Michigan)/CareQOL Project/"
}
confint.geeglm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(cc),
                cbind(est=round(Estimate,3),
                      lwr=Estimate-mult*Std.err,
                      upr=Estimate+mult*Std.err,
                      pvalue=`Pr(>|W|)`))
  rownames(citab) <- rownames(cc)
  citab[parm,]
}
my_forest_plot = function(fit,sym=FALSE,var_int=NULL,force_limit=NULL,return_yabsmax=FALSE){
  d = confint(fit) %>% data.frame() %>% rownames_to_column("term") %>% filter(term!="(Intercept)") %>% mutate(shape = if_else(pvalue<0.05,8,16))
  yabsmax = NULL
  plt<-ggplot(data=d, 
              aes(x=term, y=est, ymin=lwr, ymax=upr,color=(est>=0.00),shape=shape)) +
    scale_shape_identity() +
    geom_pointrange(size=1) + 
    geom_hline(yintercept=0, lty=2) + 
    xlab("") + ylab("Mean (95% CI)") +
    theme_bw() +
    theme(legend.position = "none",text = element_text(size = 18)) +
    scale_x_discrete(limits=rev) + 
    ggtitle(as.list(attr(fit$terms,"predvars"))[[2]])
  if (sym){
    yabsmax = max(abs(d %>% select(lwr,upr)))
  }
  if (!is.null(var_int)){
    yabsmax = max(abs(d %>% filter(term == var_int) %>% select(lwr,upr))) * 2
  }
  if (!is.null(force_limit)){
    yabsmax = force_limit
  }
  
  if(!is.null(yabsmax)){
    plt = plt + coord_flip(ylim=c(-yabsmax,yabsmax))
  }else{
    plt = plt + coord_flip()
  }
  
  if (return_yabsmax)
    return(yabsmax)
  else
    return(plt)
}

calc_se = function(fit,input_trt,input_bl){
  X1 = as.matrix(model.matrix(delete.response(fit$terms),input_trt,xlev=fit$xlevels))
  X0 = as.matrix(model.matrix(delete.response(fit$terms),input_bl,xlev=fit$xlevels))
  ses = c()
  for (i in 1:dim(X1)[1]){
    se = sqrt(t(matrix(X1[i,] - X0)) %*% vcov(fit) %*% matrix(X1[i,] - X0))
    ses = c(ses,se)
  }
  return(ses)
}

my_reflection_plot = function(fit,trt_name,trt_levels,mod_name,log=FALSE,reflect_point=FALSE,week=FALSE,q1=0.25,q2=0.75){
  res = tibble()
  if (is.numeric(trt_levels[1]))
    bl_level = 0
  else
    bl_level = "None"
  hist_dat = fit$data %>% pull(!!mod_name)
  mod_range = seq(from=quantile(hist_dat,q1),to=quantile(hist_dat,q2),length.out=100)
  hist_dat = hist_dat[hist_dat <= mod_range[length(mod_range)] & hist_dat >= mod_range[1]]
  for (value in mod_range){
    if ( week != TRUE){
      input_trt = expand_grid(!!trt_name := trt_levels, !!mod_name := value,week=0,dem_sex="female",dem_age=60)
      input_bl = expand_grid(!!trt_name := bl_level, !!mod_name := value,week=0,dem_sex="female",dem_age=60)
    }else{
      input_trt = expand_grid(!!trt_name := trt_levels, !!mod_name := value,dem_sex="female",dem_age=60)
      input_bl = expand_grid(!!trt_name := bl_level, !!mod_name := value,dem_sex="female",dem_age=60)
    }
    if (log == TRUE){
      effects = predict(fit,input_trt)-predict(fit,input_bl)
      ses = calc_se(fit,input_trt,input_bl)
      effects_ub = effects + 1.96 * ses
      effects_lb = effects - 1.96 * ses
    }else{
      effect = predict(fit,input_trt) - predict(fit,input_bl)
      ses = calc_se(fit,input_trt,input_bl)
      effects_ub = effects + 1.96 * ses
      effects_lb = effects - 1.96 * ses
    }
    names(effects) = trt_levels
    res <- res %>% bind_rows(data.frame(effect=effects,effect_ub=effects_ub,effect_lb=effects_lb) %>% rownames_to_column("trt_level") %>% mutate(mod_value = value))
  }
  v = fit$data %>% pull(!!mod_name) %>% mean
  
  res %>% ggplot() + 
    geom_hline(yintercept = 0, linetype="dashed") + # geom_vline(xintercept = v,linetype="dashed") + 
    geom_line(aes(x=mod_value,y=effect,linetype=trt_level),size=1.0) +
    geom_ribbon(aes(x=mod_value,ymin=effect_lb,ymax=effect_ub,linetype=trt_level),alpha=0.0,color="grey",size=1.0)+#geom_freqpoly(aes(x=mod,y=..count../50*diff(range(res$effect))/4+min(res$effect)),position="identity",data=data.frame(mod = hist_dat))+
    scale_y_continuous(name="Effect",sec.axis = sec_axis(~. ,name="")) +
    scale_linetype_manual(values=c("solid","dashed","dotted")) + 
    theme_bw() + #geom_text(data=data.frame(x=v,y=min(res$effect)), aes(x, y), label=floor(v),vjust=1) +
    theme(axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank()) +
    xlab(mod_name)+
    theme(legend.position = "bottom")-> p
  ylim_axis = layer_scales(p)$y$get_limits()
  p = p + geom_segment(aes(x=mod,xend=mod,y=ylim_axis[1],yend=1/5*diff(ylim_axis)/4+ylim_axis[1]),data=data.frame(mod = hist_dat))
  if (reflect_point == TRUE){
    v = exp(-coef(fit)[trt_name]/coef(fit)[paste(trt_name,":log(",mod_name,")",sep="")])
    p = p + geom_vline(xintercept = v, linetype="dashed") + geom_text(data=data.frame(x=v,y=min(res$effect)), aes(x, y), label=floor(v),vjust=1)
  }
  # p + theme(legend.position = c(0.12,0.92),legend.background = element_rect(fill = "white", color = "black"))
  p + theme(legend.position = "none")
}

my_reflection_plot(fit=fit.stress_4lvl_step,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_step",log=TRUE)

# load data
stress_item0 = "TBICQ_S32r"
sadness_item0 = "EDDEP29"
worry_item0 = "EDANX53"
dat <- read_csv(paste(prefix,"AS_data/Full Dataset/full_data.csv",sep = ""),guess_max = 3000) %>%
  select(participant_id,date,Type,starts_with("item"),starts_with("t_score"),starts_with("t_se"),msg_sent,step_count,sleep_count) %>%
  rename(item_id_CaregiverStress="item_id_Caregiver Stress",
         t_score_CaregiverStress="t_score_Caregiver Stress",
         t_se_CaregiverStress="t_se_Caregiver Stress") %>%
  mutate(date=as.Date(date)) %>%
  mutate(weekday = weekdays(date)) %>% 
  left_join(read_csv(paste(prefix,"helper_files/participants_demo.csv",sep="")) %>% dplyr::rename(participant_id=participant_ID))
tmp = dat
week_num_assgn = data.frame(date = seq(min(dat$date),max(dat$date),1)) %>% tibble() %>%
  mutate(week = as.integer(floor((date-min(date))/7)))
tmp <- tmp %>% left_join(week_num_assgn) %>% group_by(participant_id) 
tmp1 <- tmp %>% group_by(participant_id,week) %>% tally() %>%
  filter(n==7)
tmp <- tmp1 %>% left_join(tmp) %>% mutate(day = as.integer(date-min(date)))%>% mutate(week=week-min(week))

##############################################Caregiver Stress############################################################

# inclusion criteria
filter1 = tmp %>% filter(!is.na(t_score_CaregiverStress)) %>%
  arrange(participant_id,week,date) %>%
  group_by(participant_id,week) %>%
  mutate(count_items=n()) %>%
  filter(row_number() == n())  %>%
  ungroup() %>%
  filter((count_items >= 3)) %>%
  distinct(participant_id,week) # filter based on number of items per week and final standard error of t score
anti_filter1 = tmp %>% filter(is.na(msg_sent)) %>% distinct(participant_id,week)
anti_filter2 = tmp %>% filter(is.na(dem_sex)) %>% distinct(participant_id,week)

# apply filters
tmp0 = tmp %>% anti_join(anti_filter1) %>% anti_join(anti_filter2) %>% inner_join(filter1)
tmp0_pre_measures = tmp0 %>% group_by(participant_id,week) %>% 
  summarise(pre_week_step = mean(step_count,na.rm=TRUE),
            pre_week_step_log = mean(log(step_count),na.rm=TRUE),
            pre_week_sleep = mean(sleep_count,na.rm=TRUE)) %>% 
  mutate(week = week + 1) %>% 
  ungroup() #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_pre_t_score_CaregiverStress = tmp0 %>% group_by(participant_id,week) %>% 
  filter(!is.na(t_score_CaregiverStress)) %>% 
  filter(row_number() == n()) %>%
  mutate(week = week + 1) %>% 
  ungroup() %>%
  select(participant_id,week,t_score_CaregiverStress) %>%
  rename(pre_week_t_score_CaregiverStress = t_score_CaregiverStress) %>%
  mutate(pre_week_t_score_CaregiverStress_log = log(pre_week_t_score_CaregiverStress)) #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_msg_num = tmp0 %>% mutate(final_t_score_weekday=as.integer(format(date,"%u"))) %>%
  group_by(participant_id,week) %>% 
  filter(!is.na(t_score_CaregiverStress)) %>% 
  filter(row_number() == n()) %>%
  select(participant_id,week,final_t_score_weekday) %>% 
  full_join(tmp0 %>% select(participant_id,week,msg_sent,date)%>% mutate(weekday=as.integer(format(date,"%u")))) %>%
  filter(final_t_score_weekday > weekday) %>%
  group_by(participant_id,week) %>%
  summarise(msg_num_upto_final_t_score = sum(msg_sent))


# analysis
## get the latest values given partcipant_id and week
tmp1 = tmp0 %>%
  left_join(tmp0_msg_num) %>%
  filter(!is.na(t_score_CaregiverStress)) %>% 
  group_by(participant_id,week) %>%
  filter(row_number() == n()) %>% 
  inner_join(tmp0_pre_measures) %>% 
  inner_join(tmp0_pre_t_score_CaregiverStress) %>%
  mutate(msg_high = msg_num_upto_final_t_score >=3) %>%
  mutate(msg_level = case_when(msg_num_upto_final_t_score == 0 ~ "None",
                               msg_num_upto_final_t_score <= 2 ~ "Low",
                               msg_num_upto_final_t_score <= 4 ~ "Medium",
                               msg_num_upto_final_t_score <= 6 ~ "High")) %>%
  mutate(msg_level = factor(msg_level,levels = c("None","Low","Medium","High")))
fit.stress_4lvl <- geeglm(log(t_score_CaregiverStress) ~ msg_level+week+dem_sex+dem_age,
                          data = tmp1,id = factor(participant_id),corstr = "i")

my_forest_plot(fit.stress_4lvl,var_int = "msg_levelLow")
exp(confint(fit.stress_4lvl)[c(2,3,4),c(1,2,3)])-1
msg_none_tscore = exp(predict(fit.stress_4lvl,newdata = data.frame(msg_level="None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE))))
msg_none_tscore*(exp(confint(fit.stress_4lvl)[c(2,3,4),c(1,2,3)])-1)


# previous t_score_CaregiverStress ## not siginificant but having trend
fit.stress_4lvl_stress <- geeglm(log(t_score_CaregiverStress) ~ msg_level * log(pre_week_t_score_CaregiverStress)+week+dem_sex+dem_age,
                                 data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.stress_4lvl_stress)
png("~/Downloads/stress_4lvl_stress.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.stress_4lvl_stress,trt_name="msg_level",trt_levels=c("Low","Medium","High"),mod_name = "pre_week_t_score_CaregiverStress",log=TRUE) +
  xlab("Previous week's T score of caregiver strain") +
  ylab("Effect of the JITAI")
dev.off()


# week-in-study ## neither significant nor trend
fit.stress_4lvl_week <- geeglm(log(t_score_CaregiverStress) ~ msg_level *week+dem_sex+dem_age,
                               data = tmp1 %>% mutate(week=week-1),id = factor(participant_id),corstr = "i")
summary(fit.stress_4lvl_week)
png("~/Downloads/stress_4lvl_week.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.stress_4lvl_week,trt_name="msg_level",trt_levels=c("Low","Medium","High"),mod_name = "week",log=TRUE,week=TRUE,q1=0,q2=1)+
  xlab("Week in the study") +
  ylab("Effect of the JITAI")+ 
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))
dev.off()


# previous week step count ## having trend
fit.stress_4lvl_step <- geeglm(log(t_score_CaregiverStress) ~ msg_level*log(pre_week_step)+ week+dem_sex+dem_age,
                               data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.stress_4lvl_step)
png("~/Downloads/stress_4lvl_step.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.stress_4lvl_step,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_step",log=TRUE)+
  xlab("Previous week's daily average step count") +
  ylab("Effect of the JITAI") 
dev.off()

# previous week sleep ## no trend
fit.stress_4lvl_sleep <- geeglm(log(t_score_CaregiverStress) ~ msg_level*pre_week_sleep+ week+dem_sex+dem_age,
                                data = tmp1 %>%  filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
summary(fit.stress_4lvl_sleep)
png("~/Downloads/stress_4lvl_sleep.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.stress_4lvl_sleep,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_sleep",log=TRUE)+
  xlab("Previous week's daily average sleep minutes") +
  ylab("Effect of the JITAI")
dev.off()

# main effect by different subgroups
fit.stress_4lvl_subgroup <- geeglm(log(t_score_CaregiverStress) ~ Type+msg_level:Type+week+dem_sex+dem_age,
                          data = tmp1,id = factor(participant_id),corstr = "i")
gdat <- confint(fit.stress_4lvl_subgroup)[7:15,] %>% data.frame() %>% rownames_to_column("tmp") %>%
  mutate(subgroup = case_when(str_detect(tmp,"HCT") ~ "HCT",
                              str_detect(tmp,"SCI") ~ "SCI",
                              TRUE ~ "HD"),
         msg = case_when(str_detect(tmp,"Low") ~ "Low",
                         str_detect(tmp,"Medium") ~ "Medium",
                         TRUE ~ "High")) %>%
  mutate(msg = factor(msg,levels = c("Low","Medium","High"))) %>%
  mutate(score="Caregiver strain")
gdat_subgroup = gdat
gdat %>% ggplot(aes(x=msg,y=est,color=subgroup)) + 
  theme_bw() +
  geom_point(position = position_dodge(width=0.3),size=3) + 
  geom_errorbar(aes(ymin=lwr,ymax=upr),position = position_dodge(width=0.3),width=0.2,size=1) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  xlab("Frequency of prompts") + 
  ylab("Estimate")+
  ggtitle("Caregiver Strain") + 
  scale_color_grey(start=0.0,end=0.6)

##############################################   Worry ############################################################

# inclusion criteria
filter1 = tmp %>% filter(!is.na(t_score_Worry)) %>%
  arrange(participant_id,week,date) %>%
  group_by(participant_id,week) %>%
  mutate(count_items=n()) %>%
  filter(row_number() == n())  %>%
  ungroup() %>%
  filter((count_items >= 3)) %>%
  distinct(participant_id,week) # filter based on number of items per week and final standard error of t score
anti_filter1 = tmp %>% filter(is.na(msg_sent)) %>% distinct(participant_id,week)
anti_filter2 = tmp %>% filter(is.na(dem_sex)) %>% distinct(participant_id,week)

# apply filters
tmp0 = tmp %>% anti_join(anti_filter1) %>% anti_join(anti_filter2) %>% inner_join(filter1)
tmp0_pre_measures = tmp0 %>% group_by(participant_id,week) %>% 
  summarise(pre_week_step = mean(step_count,na.rm=TRUE),
            pre_week_step_log = mean(log(step_count),na.rm=TRUE),
            pre_week_sleep = mean(sleep_count,na.rm=TRUE)) %>% 
  mutate(week = week + 1) %>% 
  ungroup() #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_pre_t_score_Worry = tmp0 %>% group_by(participant_id,week) %>% 
  filter(!is.na(t_score_Worry)) %>% 
  filter(row_number() == n()) %>%
  mutate(week = week + 1) %>% 
  ungroup() %>%
  select(participant_id,week,t_score_Worry) %>%
  rename(pre_week_t_score_Worry = t_score_Worry) %>%
  mutate(pre_week_t_score_Worry_log = log(pre_week_t_score_Worry)) #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_msg_num = tmp0 %>% mutate(final_t_score_weekday=as.integer(format(date,"%u"))) %>%
  group_by(participant_id,week) %>% 
  filter(!is.na(t_score_Worry)) %>% 
  filter(row_number() == n()) %>%
  select(participant_id,week,final_t_score_weekday) %>% 
  full_join(tmp0 %>% select(participant_id,week,msg_sent,date)%>% mutate(weekday=as.integer(format(date,"%u")))) %>%
  filter(final_t_score_weekday > weekday) %>%
  group_by(participant_id,week) %>%
  summarise(msg_num_upto_final_t_score = sum(msg_sent))

# analysis
## get the latest values given partcipant_id and week
tmp1 = tmp0 %>%
  left_join(tmp0_msg_num) %>%
  filter(!is.na(t_score_Worry)) %>% 
  group_by(participant_id,week) %>%
  filter(row_number() == n()) %>% 
  inner_join(tmp0_pre_measures) %>% 
  inner_join(tmp0_pre_t_score_Worry) %>%
  mutate(msg_high = msg_num_upto_final_t_score >=3)%>%
  mutate(msg_level = case_when(msg_num_upto_final_t_score == 0 ~ "None",
                               msg_num_upto_final_t_score <= 2 ~ "Low",
                               msg_num_upto_final_t_score <= 4 ~ "Medium",
                               msg_num_upto_final_t_score <= 6 ~ "High")) %>%
  mutate(msg_level = factor(msg_level,levels = c("None","Low","Medium","High")))
fit.worry_4lvl <- geeglm(log(t_score_Worry) ~ msg_level+ week+dem_sex+dem_age,
                         data = tmp1,id = factor(participant_id),corstr = "i")
my_forest_plot(fit.worry_4lvl)

confint(fit.worry_4lvl)[c(2,3,4),c(1,2,3)]
exp(confint(fit.worry_4lvl)[c(2,3,4),c(1,2,3)])-1
msg_none_tscore = exp(predict(fit.worry_4lvl,newdata = data.frame(msg_level="None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE))))
msg_none_tscore*(exp(confint(fit.worry_4lvl)[c(2,3,4),c(1,2,3)])-1)

# previous t_score_Worry       ## pending
fit.worry_4lvl_worry <- geeglm(log(t_score_Worry) ~ msg_level*log(pre_week_t_score_Worry) + week+dem_sex+dem_age,
                               data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.worry_4lvl_worry)
png("~/Downloads/worry_4lvl_worry.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.worry_4lvl_worry,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_t_score_Worry",log=TRUE)+
  xlab("Previous week's T score of anxiety") +
  ylab("Effect of the JITAI")
dev.off()

# week-in-study              ## pending
fit.worry_cnt_week <- geeglm(log(t_score_Worry) ~ msg_num_upto_final_t_score *week+dem_sex+dem_age,
                             data = tmp1 %>% mutate(week=week-1),id = factor(participant_id),corstr = "i")
fit.worry_4lvl_week <- geeglm(log(t_score_Worry) ~ msg_level*week+dem_sex+dem_age,
                              data = tmp1 %>% mutate(week=week-1),id = factor(participant_id),corstr = "i")
summary(fit.worry_cnt_week)
summary(fit.worry_4lvl_week)
# my_reflection_plot(fit=fit.worry_cnt_week,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "week",log=TRUE,week=TRUE)
png("~/Downloads/worry_4lvl_week.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.worry_4lvl_week,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "week",log=TRUE,week=TRUE,q1=0,q2=1)+
  xlab("Week in the study") +
  ylab("Effect of the JITAI")+ 
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))
dev.off()

# previous week step count    ## pending
fit.worry_cnt_step <- geeglm(log(t_score_Worry) ~ msg_num_upto_final_t_score *log(pre_week_step) + week+dem_sex+dem_age,
                             data = tmp1,id = factor(participant_id),corstr = "i")
fit.worry_4lvl_step <- geeglm(log(t_score_Worry) ~ msg_level*log(pre_week_step)+ week+dem_sex+dem_age,
                              data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.worry_cnt_step)
summary(fit.worry_4lvl_step)
# my_reflection_plot(fit=fit.worry_cnt_step,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "pre_week_step",log=TRUE,reflect_point = FALSE)
png("~/Downloads/worry_4lvl_step.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.worry_4lvl_step,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_step",log=TRUE)+
  xlab("Previous week's daily average step count") +
  ylab("Effect of the JITAI")
dev.off()

# previous week sleep
fit.worry_cnt_sleep <- geeglm(log(t_score_Worry) ~ msg_num_upto_final_t_score *pre_week_sleep + week+dem_sex+dem_age,
                              data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
fit.worry_4lvl_sleep <- geeglm(log(t_score_Worry) ~ msg_level*pre_week_sleep+ week+dem_sex+dem_age,
                               data = tmp1 %>%  filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
summary(fit.worry_cnt_sleep)
summary(fit.worry_4lvl_sleep)
# my_reflection_plot(fit=fit.worry_cnt_sleep,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "pre_week_sleep",log=TRUE)
png("~/Downloads/worry_4lvl_sleep.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.worry_4lvl_sleep,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_sleep",log=TRUE)+
  xlab("Previous week's daily average sleep minutes") +
  ylab("Effect of the JITAI")
dev.off()


fit.worry_4lvl_subgroup <- geeglm(log(t_score_Worry) ~ Type+msg_level:Type+week+dem_sex+dem_age,
                                   data = tmp1,id = factor(participant_id),corstr = "i")
gdat <- confint(fit.worry_4lvl_subgroup)[7:15,] %>% data.frame() %>% rownames_to_column("tmp") %>%
  mutate(subgroup = case_when(str_detect(tmp,"HCT") ~ "HCT",
                              str_detect(tmp,"SCI") ~ "SCI",
                              TRUE ~ "HD"),
         msg = case_when(str_detect(tmp,"Low") ~ "Low",
                         str_detect(tmp,"Medium") ~ "Medium",
                         TRUE ~ "High")) %>%
  mutate(msg = factor(msg,levels = c("Low","Medium","High"))) %>%
  mutate(score="Anxiety")
gdat_subgroup = rbind(gdat_subgroup,gdat)
gdat %>% ggplot(aes(x=msg,y=est,color=subgroup)) + 
  theme_bw() +
  geom_point(position = position_dodge(width=0.3),size=3) + 
  geom_errorbar(aes(ymin=lwr,ymax=upr),position = position_dodge(width=0.3),width=0.2,size=1) + 
  geom_hline(yintercept = 0, linetype="dashed")+ 
  xlab("Frequency of prompts") + 
  ylab("Estimate")+
  ggtitle("Anxiety") + 
  scale_color_grey(start=0.0,end=0.6)

##############################################Sadness############################################################

# inclusion criteria
filter1 = tmp %>% filter(!is.na(t_score_Sadness)) %>%
  arrange(participant_id,week,date) %>%
  group_by(participant_id,week) %>%
  mutate(count_items=n()) %>%
  filter(row_number() == n())  %>%
  ungroup() %>%
  filter((count_items >= 3)) %>%
  distinct(participant_id,week) # filter based on number of items per week and final standard error of t score
anti_filter1 = tmp %>% filter(is.na(msg_sent)) %>% distinct(participant_id,week)
anti_filter2 = tmp %>% filter(is.na(dem_sex)) %>% distinct(participant_id,week)

# apply filters
tmp0 = tmp %>% anti_join(anti_filter1) %>% anti_join(anti_filter2) %>% inner_join(filter1)
tmp0_pre_measures = tmp0 %>% group_by(participant_id,week) %>% 
  summarise(pre_week_step = mean(step_count,na.rm=TRUE),
            pre_week_step_log = mean(log(step_count),na.rm=TRUE),
            pre_week_sleep = mean(sleep_count,na.rm=TRUE)) %>% 
  mutate(week = week + 1) %>% 
  ungroup() #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_pre_t_score_Sadness = tmp0 %>% group_by(participant_id,week) %>% 
  filter(!is.na(t_score_Sadness)) %>% 
  filter(row_number() == n()) %>%
  mutate(week = week + 1) %>% 
  ungroup() %>%
  select(participant_id,week,t_score_Sadness) %>%
  rename(pre_week_t_score_Sadness = t_score_Sadness) %>%
  mutate(pre_week_t_score_Sadness_log = log(pre_week_t_score_Sadness)) #%>%
#mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
tmp0_msg_num = tmp0 %>% mutate(final_t_score_weekday=as.integer(format(date,"%u"))) %>%
  group_by(participant_id,week) %>% 
  filter(!is.na(t_score_Sadness)) %>% 
  filter(row_number() == n()) %>%
  select(participant_id,week,final_t_score_weekday) %>% 
  full_join(tmp0 %>% select(participant_id,week,msg_sent,date)%>% mutate(weekday=as.integer(format(date,"%u")))) %>%
  filter(final_t_score_weekday > weekday) %>%
  group_by(participant_id,week) %>%
  summarise(msg_num_upto_final_t_score = sum(msg_sent))

# analysis
## get the latest values given partcipant_id and week
tmp1 = tmp0 %>%
  left_join(tmp0_msg_num) %>%
  filter(!is.na(t_score_Sadness)) %>% 
  group_by(participant_id,week) %>%
  filter(row_number() == n()) %>% 
  inner_join(tmp0_pre_measures) %>% 
  inner_join(tmp0_pre_t_score_Sadness) %>%
  mutate(msg_high = msg_num_upto_final_t_score >=3)%>%
  mutate(msg_level = case_when(msg_num_upto_final_t_score == 0 ~ "None",
                               msg_num_upto_final_t_score <= 2 ~ "Low",
                               msg_num_upto_final_t_score <= 4 ~ "Medium",
                               msg_num_upto_final_t_score <= 6 ~ "High")) %>%
  mutate(msg_level = factor(msg_level,levels = c("None","Low","Medium","High")))

fit.sadness_4lvl <- geeglm(log(t_score_Sadness) ~ msg_level+ week+dem_sex+dem_age,
                           data = tmp1,id = factor(participant_id),corstr = "i")
my_forest_plot(fit.sadness_4lvl)
exp(confint(fit.sadness_4lvl)[c(2,3,4),c(1,2,3)])-1
msg_none_tscore = exp(predict(fit.sadness_4lvl,newdata = data.frame(msg_level="None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE))))
msg_none_tscore*(exp(confint(fit.sadness_4lvl)[c(2,3,4),c(1,2,3)])-1)

# previous t_score_Sadness
fit.sadness_cnt_sadness <- geeglm(log(t_score_Sadness) ~ msg_num_upto_final_t_score *log(pre_week_t_score_Sadness) + week+dem_sex+dem_age,
                                  data = tmp1,id = factor(participant_id),corstr = "i")
fit.sadness_4lvl_sadness <- geeglm(log(t_score_Sadness) ~ msg_level*log(pre_week_t_score_Sadness) + week+dem_sex+dem_age,
                                   data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.sadness_cnt_sadness)
summary(fit.sadness_4lvl_sadness)
png("~/Downloads/sadness_4lvl_sadness.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.sadness_4lvl_sadness,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_t_score_Sadness",log=TRUE)+
  xlab("Previous week's T score of depression") +
  ylab("Effect of the JITAI")
dev.off()
# my_reflection_plot(fit=fit.sadness_cnt_sadness,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "pre_week_t_score_Sadness",log=TRUE)
(exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "High",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=50))-predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=50)))-1) * exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=50)))
exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "High",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=40)))-exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=40)))
exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "High",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=60)))-exp(predict(fit.sadness_4lvl_sadness,newdata = data.frame(msg_level = "None",week=0,dem_sex="female",dem_age=mean(tmp1$dem_age,na.rm=TRUE),pre_week_t_score_Sadness=60)))



# week-in-study
fit.sadness_cnt_week <- geeglm(log(t_score_Sadness) ~ msg_num_upto_final_t_score *week+dem_sex+dem_age,
                               data = tmp1 %>% mutate(week=week-1),id = factor(participant_id),corstr = "i")
fit.sadness_4lvl_week <- geeglm(log(t_score_Sadness) ~ msg_level*week+dem_sex+dem_age,
                                data = tmp1 %>% mutate(week=week-1),id = factor(participant_id),corstr = "i")
summary(fit.sadness_cnt_week)
summary(fit.sadness_4lvl_week)
# my_reflection_plot(fit=fit.sadness_cnt_week,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "week",log=TRUE,week = TRUE)
png("~/Downloads/sadness_4lvl_week.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.sadness_4lvl_week,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "week",log=TRUE,week=TRUE,q1=0,q2=1)+
  xlab("Week in the study") +
  ylab("Effect of the JITAI") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))
dev.off()

# previous week step count
fit.sadness_cnt_step <- geeglm(log(t_score_Sadness) ~ msg_num_upto_final_t_score *log(pre_week_step) + week+dem_sex+dem_age,
                               data = tmp1,id = factor(participant_id),corstr = "i")
fit.sadness_4lvl_step <- geeglm(log(t_score_Sadness) ~ msg_level*log(pre_week_step)+ week+dem_sex+dem_age,
                                data = tmp1,id = factor(participant_id),corstr = "i")
summary(fit.sadness_cnt_step)
summary(fit.sadness_4lvl_step)
# my_reflection_plot(fit=fit.sadness_cnt_step,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "pre_week_step",log=TRUE,reflect_point = FALSE)
png("~/Downloads/sadness_4lvl_step.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.sadness_4lvl_step,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_step",log=TRUE)+
  xlab("Previous week's daily average step count") +
  ylab("Effect of the JITAI")
dev.off()

# previous week sleep
fit.sadness_cnt_sleep <- geeglm(log(t_score_Sadness) ~ msg_num_upto_final_t_score *pre_week_sleep + week+dem_sex+dem_age,
                                data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
fit.sadness_4lvl_sleep <- geeglm(log(t_score_Sadness) ~ msg_level*pre_week_sleep+ week+dem_sex+dem_age,
                                 data = tmp1 %>%  filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
summary(fit.sadness_cnt_sleep)
summary(fit.sadness_4lvl_sleep)
# my_reflection_plot(fit=fit.sadness_cnt_sleep,trt_name = "msg_num_upto_final_t_score",trt_levels = c(1,2,3,4,5,6),mod_name = "pre_week_sleep",log=TRUE)
png("~/Downloads/sadness_4lvl_sleep.png",width = 1000,height = 900,res=200)
my_reflection_plot(fit=fit.sadness_4lvl_sleep,trt_name = "msg_level",trt_levels = c("Low","Medium","High"),mod_name = "pre_week_sleep",log=TRUE)+
  xlab("Previous week's daily average sleep minutes") +
  ylab("Effect of the JITAI")
dev.off()

fit.sadness_4lvl_subgroup <- geeglm(log(t_score_Sadness) ~ Type+msg_level:Type+week+dem_sex+dem_age,
                                  data = tmp1,id = factor(participant_id),corstr = "i")
gdat <- confint(fit.sadness_4lvl_subgroup)[7:15,] %>% data.frame() %>% rownames_to_column("tmp") %>%
  mutate(subgroup = case_when(str_detect(tmp,"HCT") ~ "HCT",
                              str_detect(tmp,"SCI") ~ "SCI",
                              TRUE ~ "HD"),
         msg = case_when(str_detect(tmp,"Low") ~ "Low",
                         str_detect(tmp,"Medium") ~ "Medium",
                         TRUE ~ "High")) %>%
  mutate(msg = factor(msg,levels = c("Low","Medium","High")))%>%
  mutate(score="Depression")
gdat_subgroup = rbind(gdat_subgroup,gdat)
gdat %>% ggplot(aes(x=msg,y=est,color=subgroup)) + 
  theme_bw() +
  geom_point(position = position_dodge(width=0.3),size=3) + 
  geom_errorbar(aes(ymin=lwr,ymax=upr),position = position_dodge(width=0.3),width=0.2,size=1) + 
  geom_hline(yintercept = 0, linetype="dashed")+ 
  xlab("Frequency of prompts") + 
  ylab("Estimate")+
  ggtitle("Depression") + 
  scale_color_grey(start=0.0,end=0.6)


#######################
png("~/../Downloads/subgroup_plot.png",width = 2000,height = 900,res=200)
gdat_subgroup %>% ggplot(aes(x=msg,y=est,color=subgroup)) + 
  theme_bw() +
  geom_point(position = position_dodge(width=0.3),size=2) + 
  geom_errorbar(aes(ymin=lwr,ymax=upr),position = position_dodge(width=0.3),width=0.2,size=1) + 
  geom_hline(yintercept = 0, linetype="dashed")+ 
  facet_wrap(~score) + xlab("Frequency of prompts") + 
  ylab("Estimate")+
  scale_color_grey(start=0.0,end=0.7) +
  theme(text = element_text(size = 15))
dev.off()

############################################
################ QIC #######################
############################################
# for (f in list(fit.stress_cnt,fit.stress_ord,fit.stress_bin)){
#   print(QIC(f))
# }
# 
# for (f in list(fit.worry_cnt,fit.worry_ord,fit.worry_bin)){
#   print(QIC(f))
# }
# 
# for (f in list(fit.sadness_cnt,fit.sadness_ord,fit.sadness_bin)){
#   print(QIC(f))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #######################
# # supplementary plots
# tmp1 %>% ggplot(aes(y=msg_num_upto_final_t_score)) + geom_bar(stat = "count") +
#   stat_count(aes(label = ..count..),geom = "text",position=position_stack(vjust=1.06),colour = "blue") + 
#   scale_y_reverse(breaks = 0:7)
# 
# 
# # sensitivity plot
# tmp1 <- data.frame(est=c(-0.23,-0.51,-1.02),
#                    lwr=c(-0.97,-1.23,-1.78),
#                    upr=c(0.51,0.27,-0.26),
#                    t_se_threshold=c(3,4,5))
# ggplot(data=tmp1, 
#        aes(x=t_se_threshold, y=est, ymin=lwr, ymax=upr)) +
#   geom_pointrange(size=1) + 
#   geom_hline(yintercept=0, lty=2) + 
#   xlab("") + ylab("Mean (95% CI) of the effect of weekly msg num") +
#   theme_bw() +
#   theme(legend.position = "none") + 
#   scale_x_continuous(breaks=c(3,4,5))
# 
# tmp1 = tmp %>% arrange(participant_id,date) %>% 
#   filter(!is.na(t_se_CaregiverStress) & !is.na(msg_sent)) %>% 
#   group_by(participant_id,week) %>% filter(row_number() == n())
# tmp0 %>% arrange(participant_id,date) %>% 
#   filter(!is.na(t_se_CaregiverStress) & !is.na(msg_sent)) %>% 
#   group_by(participant_id,week) %>% filter(row_number() == n())
# ggplot(tmp1) + geom_point(aes(x=t_score_CaregiverStress,y=t_se_CaregiverStress)) + 
#   theme_bw() + geom_hline(yintercept = 3,linetype="dashed",color="red")
# 
# ##############################################Caregiver Stress############################################################
# 
# # inclusion criteria
# filter1 = tmp %>% filter(!is.na(t_score_CaregiverStress)) %>%
#   arrange(participant_id,week,date) %>%
#   group_by(participant_id,week) %>%
#   mutate(count_items=n()) %>%
#   filter(row_number() == n())  %>%
#   ungroup() %>%
#   filter((count_items >= 3)) %>%
#   distinct(participant_id,week) # filter based on number of items per week and final standard error of t score
# anti_filter1 = tmp %>% filter(is.na(msg_sent)) %>% distinct(participant_id,week)
# anti_filter2 = tmp %>% filter(is.na(dem_sex)) %>% distinct(participant_id,week)
# 
# # apply filters
# tmp0 = tmp %>% anti_join(anti_filter1) %>% anti_join(anti_filter2) %>% inner_join(filter1)
# tmp0_pre_measures = tmp0 %>% group_by(participant_id,week) %>% 
#   summarise(pre_week_step = mean(step_count,na.rm=TRUE),
#             pre_week_step_log = mean(log(step_count),na.rm=TRUE),
#             pre_week_sleep = mean(sleep_count,na.rm=TRUE)) %>% 
#   mutate(week = week + 1) %>% 
#   ungroup() %>%
#   mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
# tmp0_pre_t_score_CaregiverStress = tmp0 %>% group_by(participant_id,week) %>% 
#   filter(!is.na(t_score_CaregiverStress)) %>% 
#   filter(row_number() == n()) %>%
#   mutate(week = week + 1) %>% 
#   ungroup() %>%
#   select(participant_id,week,t_score_CaregiverStress) %>%
#   rename(pre_week_t_score_CaregiverStress = t_score_CaregiverStress) %>%
#   mutate(pre_week_t_score_CaregiverStress_log = log(pre_week_t_score_CaregiverStress)) %>%
#   mutate(across(starts_with("pre"),~.x-mean(.x,na.rm=TRUE)))
# tmp0_msg_num = tmp0 %>% mutate(final_t_score_weekday=as.integer(format(date,"%u"))) %>%
#   group_by(participant_id,week) %>% 
#   filter(!is.na(t_score_CaregiverStress)) %>% 
#   filter(row_number() == n()) %>%
#   select(participant_id,week,final_t_score_weekday) %>% 
#   full_join(tmp0 %>% select(participant_id,week,msg_sent,date)%>% mutate(weekday=as.integer(format(date,"%u")))) %>%
#   filter(final_t_score_weekday > weekday) %>%
#   group_by(participant_id,week) %>%
#   summarise(msg_num_upto_final_t_score = sum(msg_sent))
# 
# 
# # analysis
# ## get the latest values given partcipant_id and week
# tmp1 = tmp0 %>%
#   left_join(tmp0_msg_num) %>%
#   filter(!is.na(t_score_CaregiverStress)) %>% 
#   group_by(participant_id,week) %>%
#   filter(row_number() == n()) %>% 
#   inner_join(tmp0_pre_measures) %>% 
#   inner_join(tmp0_pre_t_score_CaregiverStress) %>%
#   mutate(msg_high = msg_num_upto_final_t_score >=3) %>% 
#   ungroup() %>%
#   mutate(week=week-mean(week))
# fit.stress_cnt_test <- geeglm(log(t_score_CaregiverStress) ~ msg_num_upto_final_t_score*week+
#                                 msg_num_upto_final_t_score*pre_week_t_score_CaregiverStress+
#                                 msg_num_upto_final_t_score*pre_week_step_log+
#                                 msg_num_upto_final_t_score*pre_week_sleep+
#                                 dem_sex+dem_age,
#                               data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# fit.stress_ord_test <- geeglm(log(t_score_CaregiverStress) ~ factor(msg_num_upto_final_t_score)*week+
#                                 factor(msg_num_upto_final_t_score)*pre_week_t_score_CaregiverStress+
#                                 factor(msg_num_upto_final_t_score)*pre_week_step_log+
#                                 factor(msg_num_upto_final_t_score)*pre_week_sleep+
#                                 dem_sex+dem_age,
#                               data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# fit.stress_bin_test <- geeglm(log(t_score_CaregiverStress) ~ msg_high*week+
#                                 msg_high*pre_week_t_score_CaregiverStress+
#                                 msg_high*pre_week_step_log+
#                                 msg_high*pre_week_sleep+
#                                 dem_sex+dem_age,
#                               data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# my_forest_plot(fit.stress_cnt_test)
# my_forest_plot(fit.stress_ord_test)
# my_forest_plot(fit.stress_bin_test)
# ggarrange(my_forest_plot(fit.stress_cnt_test),my_forest_plot(fit.stress_ord_test),my_forest_plot(fit.stress_bin_test))
# 
# 
# fit.stress_cnt_test2 <- geeglm(log(t_score_CaregiverStress) ~ msg_num_upto_final_t_score+week+pre_week_t_score_CaregiverStress+pre_week_step_log+pre_week_sleep+dem_sex+dem_age,
#                                data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# fit.stress_ord_test2 <- geeglm(log(t_score_CaregiverStress) ~ factor(msg_num_upto_final_t_score)+week+pre_week_t_score_CaregiverStress+pre_week_step_log+pre_week_sleep+
#                                  dem_sex+dem_age,
#                                data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# fit.stress_bin_test2 <- geeglm(log(t_score_CaregiverStress) ~ msg_high+week+pre_week_t_score_CaregiverStress+pre_week_step_log+pre_week_sleep+
#                                  dem_sex+dem_age,
#                                data = tmp1 %>% filter(!is.na(pre_week_sleep)),id = factor(participant_id),corstr = "i")
# my_forest_plot(fit.stress_cnt_test2)
# my_forest_plot(fit.stress_ord_test2)
# my_forest_plot(fit.stress_bin_test2)
# ggarrange(my_forest_plot(fit.stress_cnt_test2),my_forest_plot(fit.stress_ord_test2),my_forest_plot(fit.stress_bin_test2))
# 
# 
