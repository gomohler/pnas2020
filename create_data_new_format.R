library(lubridate)
library(tidyverse)
#setwd("current directory")


dts=ymd(seq(mdy("01-22-2020"),mdy("03-31-2020"),1))
dts=format(as.Date(dts, '%Y-%m-%d'), "%m-%d-%Y")

i=1
data_state=read_csv(url(paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/",dts[i],".csv")))
data_state=data_state[,c("Province/State","Country/Region","Deaths","Confirmed","Recovered")]
data_state$date=dts[i]
for(i in 2:length(dts)){
  data_tmp=read_csv(url(paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/",dts[i],".csv")))
  if(i<61){
    data_tmp=data_tmp[,c("Province/State","Country/Region","Deaths","Confirmed","Recovered")]
    data_tmp$date=dts[i]
  }else{
    data_tmp=data_tmp[,c("Province_State","Country_Region","Deaths","Confirmed","Recovered")]
    data_tmp$date=dts[i]
    names(data_tmp)=names(data_state)
  }
  
  
  data_state=rbind(data_state,data_tmp)
  
}

state_list=c("New York","California","Indiana")

names(data_state)[1]="Province.State"
names(data_state)[2]="Country.Region"
data_state$Deaths[is.na(data_state$Deaths)]=0
data_state$Confirmed[is.na(data_state$Confirmed)]=0
agg_state=aggregate(Deaths~Province.State+date,data=data_state,FUN=sum)
output_state=spread(agg_state,date,Deaths)
output_state=output_state[is.element(output_state$Province.State,state_list),]
output_state[is.na(output_state)]=0


N=ncol(output_state)
output_tmp=output_state
output_state[,2:(N-1)]=output_state[,3:N]-output_state[,2:(N-1)]
output_state=output_state[,1:(N-1)]

write.csv(output_state[,1:ncol(output_state)],"data_branching_process/state_covid_mortality.csv",row.names=F)
write.csv(output_state[,2:ncol(output_state)],"data_branching_process/state_matlab_covid_mortality.csv",row.names=F)



agg_state=aggregate(Confirmed~Province.State+date,data=data_state,FUN=sum)
output_state=spread(agg_state,date,Confirmed)
output_state=output_state[is.element(output_state$Province.State,state_list),]
output_state[is.na(output_state)]=0


N=ncol(output_state)
output_tmp=output_state
output_state[,2:(N-1)]=output_state[,3:N]-output_state[,2:(N-1)]
output_state=output_state[,1:(N-1)]

output_state[output_state<0]=0
write.csv(output_state[,1:ncol(output_state)],"data_branching_process/state_covid_confirmed.csv",row.names=F)
write.csv(output_state[,2:ncol(output_state)],"data_branching_process/state_matlab_covid_confirmed.csv",row.names=F)


