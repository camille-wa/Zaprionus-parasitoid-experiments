library(broom)
library(data.table)
library(ggplot2)
library(lubridate)
library(ggsignif)
library(cowplot)
library(lme4)
library(emmeans)
library(gt)
library(ggpubfigs)
library(tidyr)
library(dplyr)
library(ggpubfigs)

##########################################
# Z. indianus / D. hydei competition counts 
##########################################

#read data  
ZH <- read.csv("ZHcomp.csv")
#convert to data table 
ZHcomp <- as.data.table(ZH)
# melt data 
ZHcomp.melt <- melt(ZHcomp, id.vars=c("vial","start.date","start.num","type","wasp"),measure.vars=names(ZHcomp)[6:ncol(ZHcomp)])
# split columns 
ZHcomp.melt[,end.date:=tstrsplit(variable,split="_")[[1]]]
ZHcomp.melt[,species:=tstrsplit(variable,split="_")[[2]]]
# make readable dates
ZHcomp.melt[,start.date:=mdy(start.date)]
ZHcomp.melt[,end.date:=mdy(variable)]
# time elapsed
ZHcomp.melt[,time:=end.date-start.date]
# total emerged
ZHcomp.total <- ZHcomp.melt[,.(total_emerged=sum(value,na.rm=T)),.(type,vial,start.num,wasp,species)]
# total of each species 
ZHcomp.total[type %in% c("ZH"), actstartnum := start.num / 2]
ZHcomp.total[type %in% c("Z"), actstartnum := start.num]
# proportion alive 
ZHcomp.total[,prop:=total_emerged/actstartnum]
# interspecific subset 
subset_inter <- ZHcomp.total[ZHcomp.total$type == "ZH", ]
#low and high subsets for contrasts
lowZHinter_subset <- subset_inter[subset_inter$start.num == 50,]
highZHinter_subset <- subset_inter[subset_inter$start.num == 150,]

##########################################
# Z. indianus / D. simulans competition counts 
##########################################

#read in data 
ZS <- read.csv("ZScomp_2024.csv")
#convert to data table 
ZScomp <- as.data.table(ZS)
#melt data 
ZScomp.melt <- melt(ZScomp, id.vars=c("vial","start.date","start.num","type","wasp"),measure.vars=names(ZScomp)[6:ncol(ZScomp)])
#split columns 
ZScomp.melt[,end.date:=tstrsplit(variable,split="_")[[1]]]
ZScomp.melt[,species:=tstrsplit(variable,split="_")[[2]]]
#readable dates 
ZScomp.melt[,start.date:=mdy(start.date)]
ZScomp.melt[,end.date:=mdy(variable)]
#time elapsed
ZScomp.melt[,time:=end.date-start.date]
#total emerged
ZScomp.total <- ZScomp.melt[,.(total_emerged=sum(value,na.rm=T)),.(start.date,type,vial,start.num,wasp,species)]
#total of each species(ZH)
ZScomp.total[type %in% c("ZS"), actstartnum := start.num / 2]
ZScomp.total[type %in% c("Z"), actstartnum := start.num]
#prop alive 
ZScomp.total[,prop:=total_emerged/actstartnum]
#inter subset
subset_inter_s <- ZScomp.total[ZScomp.total$type == "ZS", ]
#make low and high subsets for contrasts 
lowZSinter <- subset_inter_s[subset_inter_s$start.num == 50,]
highZSinter <- subset_inter_s[subset_inter_s$start.num == 150,]

##########################################
# Z. indianus intraspecific competition counts
##########################################

#make intra subset
Zintra <- ZHcomp.total[ZHcomp.total$species == "Z"]
subset_intra <- Zintra[Zintra$type == "Z" ]

Zintra_s <- ZScomp.total[ZScomp.total$species == "Z"]
subset_intra_s <- Zintra_s[Zintra_s$type == "Z" ]

############## combining low density ZH/ZS intra
ZHlow <- subset(subset_intra, start.num == 50)
ZSlow <- subset(Zintra_s, type == "Z" & start.num == 50)
combined_low <- rbind(ZSlow, ZHlow, fill=TRUE)

#combining high density
ZHhigh <- subset(subset_intra, start.num == 150)
ZShigh <- subset(Zintra_s, type == "Z" & start.num == 150)

combined_high <- rbind(ZShigh, ZHhigh, fill=TRUE)

#add density columns 
combined_high$density <- "high"
combined_low$density <- "low"

#combine 
intraspecific_combined <- rbind(combined_high, combined_low, fill=TRUE)

intraspecific_combined[,density:=factor(density,levels=c("low","high"))]
#make density a factor
intraspecific_combined$density <- as.factor(intraspecific_combined$density)

#######################
# competition statistical analyses 
#######################

# Intraspecific GLM
combinedintra.model <- glm(cbind(total_emerged, start.num - total_emerged) ~ (wasp*density), 
                           family = binomial, 
                           data = intraspecific_combined)
summary(combinedintra.model)

# Intraspecific contrasts
emm_combinedintra <- emmeans(combinedintra.model, ~ wasp*density)
intra.contrast <- contrast(emm_combinedintra, method = "pairwise")
print(intra.contrast)

## Z. indianus x D. hydei 
#low density GLM ZH 
lowZHinter_subset <- subset_inter[subset_inter$start.num == 50,]
lowZHinter.model <- glm(cbind(total_emerged, actstartnum - total_emerged) ~ (wasp*species), 
                        family = quasibinomial, 
                        data = lowZHinter_subset)
summary(lowZHinter.model)

#low density contrasts ZH 
emm_lowZHinter <- emmeans(lowZHinter.model, ~ species*wasp)
lowZH.contrast <- contrast(emm_lowZHinter, method = "pairwise")
print(lowZH.contrast)

#high density GLM ZH
highZHinter_subset <- subset_inter[subset_inter$start.num == 150,]

highZHinter.model <- glm(cbind(total_emerged, actstartnum - total_emerged) ~ (wasp*species), 
                         family = quasibinomial, 
                         data = highZHinter_subset)
summary(highZHinter.model)
#high density contrasts ZH 
emm_highZHinter <- emmeans(ZHinter3.model, ~ species*wasp)
highZH.contrast <- contrast(emm_highZHinter, method = "pairwise")
print(highZH.contrast) 

##### Z. indianus x D. simulans  

##low density GLM ZS
lowZSinter <- subset_inter_s[subset_inter_s$start.num == 50,]
lowZSinter.model <- glm(cbind(total_emerged, actstartnum - total_emerged) ~ (wasp*species), 
                        family = quasibinomial, 
                        data = lowZSinter)
summary(lowZSinter.model)
#low density contrasts ZS
emm_lowZSinter <- emmeans(lowZSinter.model, ~ species*wasp)
lowZS.contrast <- contrast(emm_lowZSinter, method = "pairwise")
print(lowZS.contrast)

#high density GLM ZS
highZSinter <- subset_inter_s[subset_inter_s$start.num == 150,]

highZSinter.model <- glm(cbind(total_emerged, actstartnum - total_emerged) ~ (wasp*species), 
                         family = quasibinomial, 
                         data = highZSinter)
summary(highZSinter.model)
#high density contrasts ZS 
emm_highZSinter <- emmeans(highZSinter.model, ~ species*wasp)
highZS.contrast <- contrast(emm_highZSinter, method = "pairwise")
print(highZS.contrast) 

#######################
# plotting 
#######################

#Z. indianus x D. hydei plot 

ZHinterplot <- ggplot(data = subset_inter) + 
  geom_boxplot(aes(x = wasp, y = prop, fill = species), size = 0.3) +
  facet_wrap(~ start.num, ncol = 1) +  # Stack facets vertically
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5, 7)], 
                    labels = c("Z" = expression(italic("Z. indianus")),
                               "H" = expression(italic("D. hydei")))) +  
  scale_x_discrete(labels = c("N" = "Absent", "Y" = "Present")) + 
  theme(legend.text = element_text(face = "italic"),
        panel.background = element_rect(fill = "white", color = "gray"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "gray", fill = NA),
        legend.key = element_blank()) +  
  labs(x = "Parasitoid", y = "") +
  ylim(0, 1) 

print(ZHinterplot)

#Z. indianus x D. simulans plot 

ZSinterplot <- ggplot(data = subset_inter_s) + 
  geom_boxplot(aes(x = wasp, y = prop, fill = species), size = 0.3) +
  facet_wrap(~ start.num, ncol = 1) +  # Stack facets vertically
  scale_fill_manual(values = friendly_pal("bright_seven")[c(6, 7)], 
                    labels = c("Z" = expression(italic("Z. indianus")),
                               "S" = expression(italic("D. simulans")))) +
  scale_x_discrete(labels = c("N" = "Absent", "Y" = "Present")) + 
  theme(legend.text = element_text(face = "italic"),
        panel.background = element_rect(fill = "white", color = "gray"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "gray", fill = NA),
        legend.key = element_blank()) +  
  labs(x = "", y = "") +
  ylim(0, 1)  

print(ZSinterplot)

#Z. indianus only 

intraplot_test <- ggplot(data = intraspecific_combined) + 
  geom_boxplot(aes(x = wasp, y = prop, fill = species), 
               color = "black", width = 0.4, size = 0.3) +  # Decrease width to make the boxes narrower
  scale_fill_manual(values = friendly_pal("bright_seven")[7],  
                    labels = c("Z" = bquote(atop(italic("Z. indianus"), "(intraspecific)")))) +
  labs(x = "", y = "Proportion emerged", fill = "species") +  # Set y-axis label
  scale_x_discrete(labels = c("N" = "Absent", "Y" = "Present")) +  # Update x-axis labels
  theme(panel.background = element_rect(fill = "white", color = "gray"),
        panel.grid.major = element_blank(),
        legend.key = element_rect(color = NA)) +
  ylim(0, 1) +  # Keep the same y-axis limits
  facet_wrap(~ density, labeller = as_labeller(c("low" = "50", "high" = "150")), ncol = 1)  # Stack facets vertically

print(intraplot_test)


#Combine ZH/ZS/ZAP only interplots
ZHZSZZinterplot <- plot_grid(
  ZHinterplot, ZSinterplot, intraplot_test,
  labels = c("A", "B", "C"),    # Labels for the panels
  nrow = 1,                     # One row for horizontal layout
  align = 'h',                  # Align plots horizontally
  rel_widths = c(1, 1, 1)      # Set equal widths for each plot
)

print(ZHZSZZinterplot)

