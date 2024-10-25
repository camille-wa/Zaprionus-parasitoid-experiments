library(broom)
library(data.table)
library(ggplot2)
library(lubridate)
library(lme4)
library(dplyr)
library(ggpattern)
library(ggpubfigs)

#FL/CT 
FLCT_N <- fread ("~/Library/CloudStorage/Box-Box/FL_CT_N.csv")
FLCT_Y <- fread ("~/Library/CloudStorage/Box-Box/FL_CT_Y.csv")

#melt
FLCT_N.melt <- melt(FLCT_N, id.vars=c("vial","line", "start.date","start.num"),measure.vars=names(FLCT_N)[5:ncol(FLCT_N)])
FLCT_Y.melt <- melt(FLCT_Y, id.vars=c("vial","line", "start.date","start.num"),measure.vars=names(FLCT_Y)[5:ncol(FLCT_Y)])
#split columns 
FLCT_N.melt[,end.date:=tstrsplit(variable,split="_")[[1]]]
FLCT_Y.melt[,end.date:=tstrsplit(variable,split="_")[[1]]]
#readable dates
FLCT_N.melt[,start.date:=mdy(start.date)]
FLCT_N.melt[,end.date:=mdy(variable)]

FLCT_Y.melt[,start.date:=mdy(start.date)]
FLCT_Y.melt[,end.date:=mdy(variable)]

#time elapsed
FLCT_N.melt[,time:=end.date-start.date]
FLCT_Y.melt[,time:=end.date-start.date]

#total emerged
FLCT_N.total<- FLCT_N.melt[,.(total_emerged=sum(value,na.rm=T)),.(start.date,line,vial,start.num)]
FLCT_Y.total <- FLCT_Y.melt[,.(total_emerged=sum(value,na.rm=T)),.(start.date,line,vial,start.num)]

#prop alive 
FLCT_N.total[,prop:=total_emerged/start.num]
FLCT_Y.total[,prop:=total_emerged/start.num]

FLCT_Y.total[,state:=tstrsplit(line,split="-")[[1]]]
FLCT_Y.total[,lines:=tstrsplit(line,split="-")[[2]]]

#filter out lines 
FL_Y <- FLCT_Y.total[line %like% "FL"]
FL_N <- FLCT_N.total[line %like% "FL"]

# Subset the data table for CT lines
CT_Y<- FLCT_Y.total[line %like% "CT"]
CT_N<- FLCT_N.total[line %like% "CT"]

#GLMM 
FLCT.model <- glmer(
  cbind(total_emerged, start.num - total_emerged) ~ state + (1 | line), 
  data = FLCT_Y.total, 
  family = binomial
)

summary(FLCT.model)

#Figure 2 genetic variation plot 
ggplot(FLCT_Y.total, aes(x = state, y = prop, fill = state)) +
  geom_boxplot(fill = "white") +  
  geom_point(aes(color = line), size = 3, position = position_jitter(width = 0.2, height = 0)) +
  labs(x = "State", y = "Proportion emerged") +
  theme_minimal() +
  guides(fill = FALSE, color = FALSE)






