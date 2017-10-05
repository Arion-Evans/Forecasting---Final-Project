library(ggplot2)
library(googleVis)
library(dplyr)
library(Hmisc)
library(tidyr)
library(forcats)
library(cowplot)
library(GGally)

overall.info.1 <- overall.info
overall.info.2 <- overall.info
overall.info.1$`Model Fitted` <- fct_lump(overall.info$`Model Fitted`, n = 8)
overall.info.1$`Season Frequency` <- overall.info.1$`Season Frequency` %>% 
  factor(levels = c(1,4,12), labels = c("Yearly","Quarterly","Monthly"))
overall.info.1$`Series Category` <- ifelse(grepl("DEMO",overall.info.1$`Series Category`), "DEMOGRAPHIC", 
                                           as.character(overall.info.1$`Series Category`))
overall.info.2$`Season Frequency` <- overall.info.2$`Season Frequency` %>% 
  factor(levels = c(1,4,12), labels = c("Yearly","Quarterly","Monthly"))


# sankey
sk1 <- as.data.frame(table(overall.info.1$`Season Frequency`,overall.info.1$`Series Used`))
sk2 <- as.data.frame(table(overall.info.1$`Series Used`,overall.info.1$`Series Category`))
sk3 <- as.data.frame(table(overall.info.1$`Series Category`,overall.info.1$`Model Fitted`))

sk <- rbind(sk1,sk2,sk3)

sankey <- gvisSankey(sk, from='Var1', to='Var2', weight='Freq',
                options=list(height=550, width=1200,  title="Diagram Title"))

plot(sankey)


# bar charts

bar1 <- ggplot(data = overall.info.1, aes(x = `Season Frequency`)) + 
  geom_bar(aes(fill = `Series Used`), position = "dodge") +
  facet_wrap(~as.factor(`Series Category`))
bar1

bar2 <- ggplot(data = overall.info.1, aes(x = `Season Frequency`)) + 
  geom_bar(aes(fill = `Series Used`), position = "dodge") +
  labs(y = "Count", x = "Seasonality", title = "Series Used for Model Fitting by Seasonality") +
  theme_minimal()
bar2


seas.fit <- table(overall.info.2$`Season Frequency`,overall.info.2$`Model Fitted`) %>% prop.table(margin = 1)
seas.fit <- seas.fit*100
seas.fit <- as.data.frame(seas.fit)
seas.fit <- seas.fit[seas.fit$Freq != 0,]
bar3 <- ggplot(data = seas.fit, aes(x = Var1, y = Freq)) + 
  geom_bar(aes(fill = Var2), position = "dodge", stat = "identity") +
  labs(y = "Percent", x = "Seasonality", fill = "Model Fitted", title = "Model Fitted by Seasonality") + 
  theme_minimal() + 
  theme(legend.justification=c(1,1),legend.position=c(1,1), 
        legend.background = element_rect(fill="antiquewhite", size=0.5, linetype="solid", colour ="darkblue"),
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11))
bar3

# histograms

h1 <- ggplot(data = overall.info.monthly, aes(x = `Model MASE`)) + 
  geom_density(alpha = 1/2, fill = "mediumseagreen") +
  geom_histogram(colour = "white", bins = 30, aes(`Model MASE`,..density..), alpha = 1/2, fill = "darkslateblue") +
  geom_vline(xintercept = mean(overall.info.monthly$`Model MASE`), colour = "red", linetype = 2) + 
  annotate("text", label = "Mean", x = mean(overall.info.monthly$`Model MASE`) + 0.02, y = 3)+ 
  geom_vline(xintercept = median(overall.info.monthly$`Model MASE`), colour = "red") + 
  annotate("text", label = "Median", x = median(overall.info.monthly$`Model MASE`) - 0.02, y = 3) +
  labs(y = "Density", x = "Monthly MASE") +
  scale_x_continuous(limits = c(0,1.2))


h2 <- ggplot(data = overall.info.quarterly, aes(x = `Model MASE`)) + 
  geom_density(alpha = 1/2, fill = "mediumseagreen") +
  geom_histogram(colour = "white", bins = 30, aes(`Model MASE`,..density..), alpha = 1/2, fill = "darkslateblue") +
  geom_vline(xintercept = mean(overall.info.quarterly$`Model MASE`), colour = "red", linetype = 2) + 
  annotate("text", label = "Mean", x = mean(overall.info.quarterly$`Model MASE`) + 0.02, y = 3)+ 
  geom_vline(xintercept = median(overall.info.quarterly$`Model MASE`), colour = "red") + 
  annotate("text", label = "Median", x = median(overall.info.quarterly$`Model MASE`) - 0.02, y = 3.3) +
  labs(y = "Density", x = "Quarterly  MASE") +
  scale_x_continuous(limits = c(0,1.2))


h3 <- ggplot(data = overall.info.yearly, aes(x = `Model MASE`)) + 
  geom_density(alpha = 1/2, fill = "mediumseagreen") +
  geom_histogram(colour = "white", bins = 30, aes(`Model MASE`,..density..), alpha = 1/2, fill = "darkslateblue") +
  geom_vline(xintercept = mean(overall.info.yearly$`Model MASE`), colour = "red", linetype = 2) + 
  annotate("text", label = "Mean", x = mean(overall.info.yearly$`Model MASE`) - 0.02, y = 3)+ 
  geom_vline(xintercept = median(overall.info.yearly$`Model MASE`), colour = "red") + 
  annotate("text", label = "Median", x = median(overall.info.yearly$`Model MASE`) + 0.02, y = 4.2) +
  labs(y = "Density", x = "Yearly MASE") +
  scale_x_continuous(limits = c(0,1.2))


cow <- plot_grid(h1, h2,h3, 
          nrow = 3, ncol = 1,
          labels = c("Model MASE Distributions by Seasonality"),
          vjust = 1, hjust = -1.8, label_fontface = "bold.italic")
cow


hist2 <- ggplot(data = overall.info.1, aes(x = `Model MASE`)) + 
  geom_density(alpha = 1/2) +
  geom_histogram(colour = "white", bins = 15, aes(`Model MASE`,..density..), alpha = 1/2) +
  facet_grid(~`Series Used`)
hist2

hist3 <- ggplot(data = overall.info.1, aes(x = `Model MASE`)) + 
  geom_density(alpha = 1/2) +
  geom_histogram(colour = "white", bins = 15, aes(`Model MASE`,..density..), alpha = 1/2) +
  facet_grid(~`Series Category`)
hist3


# box plots

box1 <- ggplot(data = overall.info.1, aes(y = `Model MASE`)) + 
  geom_boxplot(aes(x = `Series Category`, fill = `Season Frequency`)) + 
  labs(fill = "Seasonality", title = "Model MASE by Series Category & Seasonality") +
  theme_minimal()
box1


# Scatter
observed.1step <- vector()
frc.1step <- vector()
seas.1step <- vector()

for(i in 1:101){
  observed.1step[i] <- series05.monthly[[i]][1]
  frc.1step[i] <- forecasts.monthly[[i]][1]
  seas.1step[i] <- "Monthly"
}
for(i in 102:202){
  observed.1step[i] <- series05.quarterly[[i-101]][1]
  frc.1step[i] <- forecasts.quarterly[[i-101]][1]
  seas.1step[i] <- "Quarterly"
}
for(i in 203:303){
  observed.1step[i] <- series05.yearly[[i-202]][1]
  frc.1step[i] <- forecasts.yearly[[i-202]][1]
  seas.1step[i] <- "Yearly"
}

df.1step <- data.frame(Observed = observed.1step,Forecasts = frc.1step, Seasonality = seas.1step)

sc <- ggplot(data = df.1step, aes(x = Observed, y = Forecasts)) + 
  stat_smooth(method = "lm",aes(colour = Seasonality)) +
  geom_point() + facet_grid(Seasonality~.) + theme_bw() + 
  guides(colour = FALSE) +
  labs(title = "1 Step Ahead Forecasts and Actual Observed Series Values by Seasonality")

sc


