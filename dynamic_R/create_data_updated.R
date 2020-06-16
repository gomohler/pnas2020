setwd("~/Dropbox/coronavirus_hawkes")

library(tidyr)

library (readr)

urlfile=urlfile="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
data<-read_csv(url(urlfile))
names(data)[1]="Province.State"
names(data)[2]="Country.Region"

N=ncol(data)

data[,5:(N-1)]=data[,6:N]-data[,5:(N-1)]
data=data[,1:N-1]

data_long <- gather(data, date,count,'1/22/20':names(data)[N-1])
data_long$date=mdy(data_long$date)



agg_data=aggregate(count~Country.Region,data=data_long,FUN=sum)
agg_data=agg_data[order(-agg_data$count),]


data_long_agg=aggregate(count~Country.Region+date,data=data_long,FUN=sum)


output=spread(data_long_agg, date, count)

clist=c("Italy","US","China")
output=output[is.element(output$Country.Region,clist),]
output=output[,1:68]

write.csv(output[,1:ncol(output)],"covid_mortality.csv",row.names=F)
write.csv(output[,2:ncol(output)],"matlab_covid_mortality.csv",row.names=F)


