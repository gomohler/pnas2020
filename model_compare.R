library(lubridate)
library(ggplot2)
library(tidyverse)
library(dplyr)
#setwd("current directory")
source('get_jhu_data.R')

source('sir.R')
source('seir.R')



st_short=c("CA","IN","NY")

output_state[output_state<0]=0
output_state_death[output_state_death<0]=0

mod_comp=data.frame(state_name=character(0),
                    data_type=character(0),
                    model=character(0),
                    rbest=numeric(0),
                    gbest=numeric(0),
                    ibest=numeric(0),
                    mubest=character(0),
                    fbest=numeric(0),
                    peak_date=character(0),
                    logL=numeric(0),
                    AIC=numeric(0))

mod_comp2=data.frame(state_name=character(0),
                    data_type=character(0),
                    model=character(0),
                    peak_date=character(0),
                    logL=numeric(0))

data_types=c("confirmed","mortality")

y=as.numeric(output_state_death[1,2:ncol(output_state_death)])  
M=length(y)
dates=seq(mdy(names(output_state)[2]),(mdy(names(output_state)[2])+200),1)
data_var_R=data.frame(date=dates[2:length(dates)])


for(k in 1:3){
  for(c_or_m in 1:2){

if(c_or_m==1){
y=as.numeric(output_state[k,2:ncol(output_state)])
}else{
y=as.numeric(output_state_death[k,2:ncol(output_state_death)])  
}  
M=length(y)
state_name=output_state$Province.State[k]
data_type=data_types[c_or_m]
dates=seq(mdy(names(output_state)[2]),(mdy(names(output_state)[2])+200),1)

rs=seq(1.5,5,.1)
gammas=seq(.01,.2,.01)

i_start=c(.005,.01,.025,.05,.1)

mus=seq(.02,.4,.02)
fs=seq(.05,1,.05)



T=200
if(k==3){
  N=19.54*1000000}
if(k==1){
  N=39.56*1000000
}
if(k==2){
  N=6.692*1000000
}

mort=.01 # mortality

best=-300000000
for(r0 in rs){
  for(is in i_start){
    for(gamma in gammas){
      
      output_sir=sir(N,T,is,0,r0,gamma,mort)
      dD=output_sir$dD
      dI=output_sir$dI
      I=output_sir$I
      
      if(c_or_m==1){
      sirPred=dI
      }else{
      sirPred=dD  
      }
      
      tbest=sum(y*log(sirPred[1:M]+10^-10))-sum(sirPred[1:M])
      RMSE=mean((sirPred[1:M]-y)^2)^.5
      
      if(r0==3&c_or_m==1){
        data_var_R$new=I
        names(data_var_R)[ncol(data_var_R)]=paste0(st_short," R=",3.0)
      }
      if(r0==2.5&c_or_m==1){
        data_var_R$new=I
        names(data_var_R)[ncol(data_var_R)]=paste0(st_short," R=",2.5)
      }
      if(r0==2&c_or_m==1){
        data_var_R$new=I
        names(data_var_R)[ncol(data_var_R)]=paste0(st_short," R=",2.0)
      }
      if(r0==1.5&c_or_m==1){
        data_var_R$new=I
        names(data_var_R)[ncol(data_var_R)]=paste0(st_short," R=",1.5)
      }
      
      ix=which.max(dI)
      peak_d=dates[ix]
      m2tmp=data.frame(state_name=state_name,
                       data_type=data_type,
                       model="sir",
                       peak_date=peak_d,
                       logL=tbest)
      
      if(abs(tbest-best)<2){
      mod_comp2=rbind(mod_comp2,m2tmp)
      }
  
      
      if(tbest>best){
        r_best=r0
        g_best=gamma
        i_best=is
        best=tbest
        dDbest=dD
        dIbest=dI
        sirPredbest=sirPred
        Ibest=I
        ix=which.max(dIbest)
        peak_best=dates[ix]
        print(c(r_best,g_best,i_best,best,peak_best,state_name))
        AIC=2*3-2*tbest
        
        mod_comp_tmp=data.frame(state_name=state_name,
                                data_type=data_type,
                            model="sir",
                            rbest=r_best,
                            gbest=g_best,
                            ibest=i_best,
                            mubest="--",
                            fbest="--",
                            peak_date=peak_best,
                            logL=tbest,
                            AIC=AIC)
        
      }
    }
  }
  }

mod_comp=rbind(mod_comp,mod_comp_tmp)

data_plot$new=sirPredbest[1:length(y)]
names(data_plot)[ncol(data_plot)]=paste0("sir_",st_short[k],"_",data_type)

best=-300000000
for(r0 in rs){
  for(is in i_start){
    for(gamma in gammas){
      for(mu in mus){
      
      output_seir=seir(N,T,is,0,r0,gamma,mu,mort)
      dD=output_seir$dD
      dI=output_seir$dI
      I=output_seir$I
      dmuE=output_seir$muE
      
      if(c_or_m==1){
        seirPred=dmuE
      }else{
        seirPred=dD  
      }
      
      tbest=sum(y*log(seirPred[1:M]+10^-10))-sum(seirPred[1:M])
      RMSE=mean((sirPred[1:M]-y)^2)^.5
      
      ix=which.max(dI)
      peak_d=dates[ix]
      m2tmp=data.frame(state_name=state_name,
                       data_type=data_type,
                       model="seir",
                       peak_date=peak_d,
                       logL=tbest)
      
      if(abs(tbest-best)<2){
        mod_comp2=rbind(mod_comp2,m2tmp)
      }
      
      if(tbest>best){
        r_best=r0
        g_best=gamma
        i_best=is
        best=tbest
        dDbest=dD
        dIbest=dI
        dmuEbest=dmuE
        Ibest=I
        mu_best=mu
        seirPredbest=seirPred
        ix=which.max(dIbest)
        peak_best=dates[ix]
        print(c(r_best,g_best,i_best,best,peak_best,state_name,mu_best))
        AIC=2*4-2*tbest
        
        mod_comp_tmp=data.frame(state_name=state_name,
                                data_type=data_type,
                                model="seir",
                                rbest=r_best,
                                gbest=g_best,
                                ibest=i_best,
                                mubest=as.character(mu_best),
                                fbest="--",
                                peak_date=peak_best,
                                logL=tbest,
                                AIC=AIC)
        
      }
      }}}
}

mod_comp=rbind(mod_comp,mod_comp_tmp)

data_plot$new=seirPredbest[1:length(y)]
names(data_plot)[ncol(data_plot)]=paste0("seir_",st_short[k],"_",data_type)
data_plot$new=y
names(data_plot)[ncol(data_plot)]=paste0("true_",st_short[k],"_",data_type)


write.csv(mod_comp,paste0("model_comparison_table",mort,".csv"),row.names=F)       

  }
}




