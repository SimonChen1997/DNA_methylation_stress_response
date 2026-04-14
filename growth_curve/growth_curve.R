library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(grid)
library(patchwork)

##########################################################################
##### import tsv
setwd("")

e_coli_od600_cfu <- read.csv("e_coli_37_1salt_LB_growth_curve.csv", header = T)
e_coli_od600_cfu$replicate <- as.character(e_coli_od600_cfu$replicate)

e_coli_37d_1nacl_LB_od600 <- read.csv("e_coli_37d_1nacl_LB_od600.csv", header = T)
e_coli_42d_1nacl_LB_od600 <- read.csv("e_coli_42d_1nacl_LB_od600.csv", header = T)
e_coli_37d_3nacl_LB_od600 <- read.csv("e_coli_37d_3nacl_LB_od600.csv", header = T)
e_coli_37d_1nacl_0p5acetic_LB_od600 <- read.csv("e_coli_37d_1nacl_0p5acetic_LB_od600.csv", header = T)
e_coli_37d_1nacl_paraffin_LB_od600 <- read.csv("e_coli_37d_1nacl_paraffin_LB_od600.csv", header = T)

##########################################################################
##### data manipulation
e_coli_od600_cfu_clean <- e_coli_od600_cfu %>% drop_na()
e_coli_od600_cfu_clean <- e_coli_od600_cfu_clean %>%
  mutate(ln_od600=log(od600))

e_coli_od600_cfu_summary <- e_coli_od600_cfu_clean %>%
  group_by(time) %>%
  summarise(mean_od600 = mean(od600),
            sd_od600 = sd(od600),
            mean_ln600 = mean(ln_od600),
            sd_ln600 = sd(ln_od600),
            mean_log10_cfu = mean(log10_cfu),
            sd_log10_cfu = sd(log10_cfu),
            mean_cfu=mean(cfu),
            sd_cfu=sd(cfu))

e_coli_od600 <- rbind(e_coli_37d_1nacl_LB_od600, e_coli_42d_1nacl_LB_od600, e_coli_37d_3nacl_LB_od600, 
                      e_coli_37d_1nacl_0p5acetic_LB_od600, e_coli_37d_1nacl_paraffin_LB_od600)

e_coli_od600_summary <- e_coli_od600 %>%
  group_by(condition, time) %>%
  summarise(mean_od600 = mean(od600),
            sd_od600 = sd(od600), .groups= "drop")

name_change_list <- list("LB_1%NaCl_37"="37", "LB_1%NaCl_42"="42", "LB_1%NaCl_37_0.5acetic"="ACE", 
                         "LB_3%NaCl_37"="TNA", "LB_1%NaCl_37_paraffin"="PA")

name_change_list_unlist <- unlist(name_change_list)
all_condition_combine$condition_short <- name_change_list_unlist[all_condition_combine$condition]
e_coli_od600_summary$condition_short <- name_change_list_unlist[e_coli_od600_summary$condition]

##########################################################################
##### function of data imputation
fit_gompertz <- function(data_table){
  pre_data_form <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("time", "pre_od600", "condition", "model", "high_time"))))
  con=c(unique(data_table$condition))
  for (c in con) {
    data= data_table[data_table$condition %in% c, ]
    i <- which.max(diff(data$mean_od600)) ## find out the highest growth rate in the data frame
    
    ## calculate the values for the formula
    a=max(data$mean_od600)
    mu=max(diff(data$mean_od600))/(data[i,"time"]-data[i-1, "time"])
    lambda=(mu*data[i, "time"]-data[i, "mean_od600"])/mu
    ## calculate the values for the formula
    starting_values <- list(a=a, ## the asymptote
                            mu=mu, ## the highest growth rate
                            lambda=lambda) ## the lag time
    print(starting_values$lambda)
    print(c)
    formula_gompertz <- "mean_od600 ~ a*exp(-exp(mu*exp(1)/a*(lambda-time)+1))"
    fitted_gompertz_model <- nlsLM(formula_gompertz, data, starting_values)
    fitted_a <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["a"]]
    fitted_mu_time <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["mu"]][["time"]]
    fitted_lambada <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["lambda"]][["time"]]
    
    pre_time=seq(0,24,by=0.5)
    y <- fitted_a*exp(-exp(fitted_mu_time*exp(1)/fitted_a*(fitted_lambada-pre_time)+1))
    mod <- "Gompertz"
    inter_dataframe <- data.frame(time=pre_time, pre_od600=y)
    high_num <- which.max(diff(inter_dataframe$pre_od600))
    high_time <- inter_dataframe[high_num, "time"]
    inter_dataframe$condition <- c
    inter_dataframe$model <- mod
    inter_dataframe$high_time <- high_time
    pre_data_form <- pre_data_form %>%
      rbind(inter_dataframe)
  }
  return(pre_data_form)
}

##########################################################################
##### data imputation
prediction_gompertz <- fit_gompertz(e_coli_od600_summary)

name_change_list <- list("LB_1%NaCl_37"="37", "LB_1%NaCl_42"="42", "LB_1%NaCl_37_0.5acetic"="ACE", 
                         "LB_3%NaCl_37"="TNA", "LB_1%NaCl_37_paraffin"="PA")

name_change_list_unlist <- unlist(name_change_list)
prediction_gompertz$condition_short <- name_change_list_unlist[prediction_gompertz$condition]
e_coli_od600_summary$condition_short <- name_change_list_unlist[e_coli_od600_summary$condition]

##########################################################################
##### od600 curve
od600 <- ggplot(e_coli_od600_cfu_summary, aes(x=time, y=mean_od600)) +
  # geom_point(color="#fd990d", size =3) +
  geom_errorbar(aes(ymin=mean_od600-sd_od600, ymax=mean_od600+sd_od600), size=0.2) +
  geom_line(color="#644b77", size=1) +
  # geom_smooth(color="#475c7a", size=0.5, se=F, formula = y ~ poly(x, 1)) +
  theme_bw() +
  labs(x="Time (h)", y="OD600") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(text = element_text(size = 20))

##### log10_cfu curve
log10_cfu <- ggplot(e_coli_od600_cfu_summary, aes(x=time, y=mean_log10_cfu)) +
  # geom_point(color= "#644b77", size =3) +
  geom_errorbar(aes(ymin=mean_log10_cfu-sd_log10_cfu, ymax=mean_log10_cfu+sd_log10_cfu), size=0.2) +
  geom_line(size=1, color="#644b77") +
  # geom_smooth(color="#475c7a", size=0.5, se=F, formula = y ~ poly(x, 1)) +
  theme_bw() +
  labs(x="Time (h)", y="Log10(CFU)") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(text = element_text(size = 20))

##### od600 - cfu
lm_od600_cfu <- ggplot(e_coli_od600_cfu_clean[e_coli_od600_cfu_clean$time > 0,], aes(x=od600, y=cfu)) +
  geom_point(color="#644b77", size =2) +
  geom_smooth(method = "lm", color="#644b77", size=1) +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), size=5, color = "#444574") +
  labs(x="OD600", y="CFU") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(text = element_text(size = 20))

predicted_normal_od600 <-ggplot() +
  geom_line(data=prediction_gompertz[prediction_gompertz$condition=="LB_1%NaCl_37",], aes(x=time, y=pre_od600, color="pre_od600"), size=1, alpha=0.9) +
  geom_line(data=e_coli_od600_cfu_summary, aes(x=time, y=mean_od600, color="mean_od600"), size=1, alpha=0.9) +
  scale_color_manual(name="Data type",
                     values=c("pre_od600"="#fd990d", "mean_od600"="#644b77"),
                     labels=c("Observed OD600", "Predicted OD600")) +
  theme_bw() +
  labs(x="Time (h)", y="OD600") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 22),
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank())

od600 <- od600 +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 22))

log10_cfu <- log10_cfu +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 22))

lm_od600_cfu <- lm_od600_cfu +
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 22))

predicted_normal_od600 <- predicted_normal_od600 +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 22))

od600 + log10_cfu + lm_od600_cfu + predicted_normal_all_con_od600 + 
  plot_layout(nrow = 2, ncol=2, heights = c(0.5,0.5), widths = c(0.5,0.5), byrow = T)

