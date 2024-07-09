library(broom)
library(data.table)
library(ggplot2)
library(lubridate)
library(lme4)
library(dplyr)
library(ggpattern)
library(ggpubfigs)

#read in data
ZSduration <- fread ("~/Library/CloudStorage/Box-Box/ZSduration.csv")

#melt data
ZSduration.melt <- melt(ZSduration, id.vars=c("vial","start.date","start.num","time","wasp"),measure.vars=names(ZSduration )[6:ncol(ZSduration )])
#split columns 
ZSduration.melt[,end.date:=tstrsplit(variable,split="_")[[1]]]
ZSduration.melt[,species:=tstrsplit(variable,split="_")[[2]]]
#readable dates
ZSduration.melt[,start.date:=mdy(start.date)]
ZSduration.melt[,end.date:=mdy(variable)]
#time elapsed
ZSduration.melt[,timeelapsed:=end.date-start.date]
#total emerged
ZSduration.total <- ZSduration.melt[,.(total_emerged=sum(value,na.rm=T)),.(start.date,time,vial,start.num,wasp,species)]
#total of each species(ZS)
ZSduration.total[, start.num := 10]
#prop alive
ZSduration.total[,prop:=total_emerged/start.num]
#filter wasp Y
ZSdurationY <- ZSduration.total[ZSduration.total$wasp == "Y", ]
ZSdurationN <- ZSduration.total[ZSduration.total$wasp == "N", ]
#add time 0 to controls
ZSdurationN <- ZSdurationN %>%
  mutate(time = 0)

#plot #add error bars 
ZSdurationN_stats <- ZSdurationN %>%
  group_by(species) %>%
  summarize(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE)
  ) 
ZSdurationY_stats <- ZSdurationY %>%
  group_by(species,time) %>%
  summarize(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE)
  )
ZSdurationY_stats <- ZSdurationY_stats %>%
  dplyr::mutate(time = as.numeric(as.character(time)))

#add time 0 to controls
ZSdurationN <- ZSdurationN %>%
  mutate(time = 0)

ZSdurationN_stats <- ZSdurationN_stats %>%
  mutate(time = 0)

ZSdurationY_stats$time <- as.factor(ZSdurationY_stats$time)
ZSdurationN_stats$time <- as.factor(ZSdurationN_stats$time)

# Combine the data frames
combined_stats <- bind_rows(ZSdurationN_stats, ZSdurationY_stats)

#combine raw data 
ZSdurationY$time <- as.factor(ZSdurationY$time)
ZSdurationN$time <- as.factor(ZSdurationN$time)

combined_raw <- bind_rows(ZSdurationN, ZSdurationY)

# add "type" column
combined_stats <- combined_stats %>%
  mutate(type = ifelse(time == 0, "control", "experimental"))

#Figure 4 host-switching plot   
ZSdurationY_plot <- ggplot(combined_stats, aes(x = as.factor(time), y = mean_prop, fill = species, pattern = type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_bar_pattern(aes(x = as.factor(time), y = mean_prop, pattern = type),
                   position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "red",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6,
                   stat = "identity") +
  geom_errorbar(aes(ymin = mean_prop, ymax = mean_prop + sd_prop),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Duration of parasitoid (Hours)",
       y = "Mean Proportion Emerged") +
  theme_minimal() +
  ggtitle("") +
  scale_fill_manual(values = friendly_pal("bright_seven")[c(6, 7)],
                    labels = c("S" = expression(italic("D. simulans")), "Z" = expression(italic("Z. indianus")))) +
  scale_pattern_manual(values = c(control = "stripe", experimental = "none")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = FALSE) + 
  ylim(0, 1)

print(ZSdurationY_plot)

#GLM 

#convert time to factor  
ZSdurationY_stats$time <- as.factor(ZSdurationY_stats$time)
ZSdurationN_stats$time <- as.factor(ZSdurationN_stats$time)
ZSdurationN_stats$time <- as.factor(as.character(ZSdurationN_stats$time))
ZSdurationY_stats$time <- as.factor(as.character(ZSdurationY_stats$time))
ZSdurationY_stats$species <- as.factor(ZSdurationY_stats$species)
ZSdurationN_stats$species <- as.factor(ZSdurationN_stats$species)
ZSdurationN_stats$species <- as.factor(as.character(ZSdurationN_stats$species))
ZSdurationY_stats$species <- as.factor(as.character(ZSdurationY_stats$species))
combined_raw$species <- as.factor(combined_raw$species)
combined_raw$time <- as.factor(combined_raw$time)

#glm 
duration.model <- glm(cbind(total_emerged, start.num - total_emerged) ~ (time*species), 
                   family = binomial, 
                   data = combined_raw)
summary(duration.model)

#emm
emm_duration <- emmeans(duration.model, ~ species * time)
#contrasts
duration.contrast <- contrast(emm_duration, method = "pairwise")
print(duration.contrast) 

emm_Z <- subset(emm_duration, species == "Z")
emm_S <- subset(emm_duration, species == "S")

z_contrasts <- contrast(emm_Z, 
                        method = list(
                          "Z 0 vs Z 1" = c(-1, 1, 0, 0),
                          "Z 0 vs Z 8" = c(-1, 0, 1, 0),
                          "Z 0 vs Z 24" = c(-1, 0, 0, 1)
                        ))

# Contrasts for species S
s_contrasts <- contrast(emm_S, 
                        method = list(
                          "S 0 vs S 1" = c(-1, 1, 0, 0),
                          "S 0 vs S 8" = c(-1, 0, 1, 0),
                          "S 0 vs S 24" = c(-1, 0, 0, 1)
                        ))

# Display results
print(z_contrasts)
print(s_contrasts)

