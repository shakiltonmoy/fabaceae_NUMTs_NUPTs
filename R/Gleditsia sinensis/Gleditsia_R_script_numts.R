setwd("D:/NUMT/R Analysis/Gleditsia sinensis")
data = read.csv("numts_assembled_columns_Gleditsia.csv")
head(data[,1:10])
i=2
datatemp = as.matrix(data[i,])
typeof(datatemp)
hist(datatemp)

startinglevel = data[i,1]
endinglevel = data[i,2]
startingcheck = matrix(1*((startinglevel - datatemp)>0))
endingcheck = matrix(1*((endinglevel - datatemp)>0))
startending_df = cbind(startingcheck,endingcheck)
dim(startending_df)
count_table = table(startending_df[,1], startending_df[,2])/1000
count_table
10*log10(count_table[2,2]/count_table[1,1])


data_nrow = nrow(data)
data_ncol = ncol(data[,-c(1,2)])
evidence = matrix(0, nrow=data_nrow, ncol=1)
for (i in 1:data_nrow){
  datatemp = as.matrix(data[i,])
  startinglevel = data[i,1]
  endinglevel = data[i,2]
  startingcheck = matrix(1*((startinglevel - datatemp)>0))
  endingcheck = matrix(1*((endinglevel - datatemp)>0))
  startending_df = cbind(startingcheck,endingcheck)
  dim(startending_df)
  count_table = table(startending_df[,1], startending_df[,2])/1000
  count_table
  evidence[i] = 10*log10(count_table[2,2]/sum(count_table[1,1],
                                              count_table[1,2],
                                              count_table[2,1]))
}
data_ev = data.frame(data, evidence=evidence)
data_ev = data_ev[order(data_ev$evidence),]
plot(data_ev$evidence)

data_ev$ev_check = 1*(data_ev$evidence>5)
table(data_ev$ev_check)

write.csv(data_ev,"Gleditsia_numts_assembled_columns_withEvidence.csv")






