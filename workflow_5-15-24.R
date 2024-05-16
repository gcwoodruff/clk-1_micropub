library(ggplot2)
library(cowplot)
library(reshape2)
library(ggforce)
library(rstatix)
library(patchwork)
library(lemon)

#set working directory
setwd("/Users/gavin//genome/clk-1")

#get data in there
dat <- read.csv("clk-1_data.csv", stringsAsFactors = TRUE)
#as factor
dat$unique_plate_id <- as.factor(dat$unique_plate_id)
#get rid of adults and gravid adults, just want total adults
datb <- subset(dat, select = -c(Adults, Gravid_Adults))


#convert wide data to long data
dat_melt <- melt(dat,measure.vars = c("Embryos","L1_L3","L4","Total_adults"))

#replace NA with 0 for plotting
dat_melt_no_NA <- dat_melt

dat_melt_no_NA[is.na(dat_melt_no_NA)] <- 0


#get inflection point and 95% confidence intervals of logistic models (i.e., time at which 50% of animals have reached a milestone)
	#define function
func_log_inf_conf <- function(dat,x,y){

	f <- function (x) 1/(1+exp(-x))
	the_glm <- glm(y ~ x, family="binomial", data=dat)
	conf <- confint(the_glm)
	#the inflection point

	p <- 0.5
	q <- (log(p/(1-p)) - coef(the_glm)[1]) / coef(the_glm)[2]

	# confidence interval at inflection point
	upper_q <- (log(p/(1-p)) - conf[1,1]) / conf[2,1]
	lower_q <- (log(p/(1-p)) - conf[1,2]) / conf[2,2]

	#return
	ret <- c(q, lower_q, upper_q)
	return(ret)
}



#again replacing NA's with 0's
datb_no_NA <- datb

datb_no_NA[is.na(datb_no_NA)] <- 0

#add L4 or older column

datb_no_NA$L4_older <- datb_no_NA$L4 + datb_no_NA$Total_adults
#add all worms per plate
datb_no_NA$total_worms_observed <-datb_no_NA$Embryos + datb_no_NA$L1_L3 + datb_no_NA$L4 + datb_no_NA$Total_adults
#fraction
datb_no_NA$fra_L1_L3 <- datb_no_NA$L1_L3/datb_no_NA$total_worms_observed
datb_no_NA$fra_L4_older <- datb_no_NA$L4_older/datb_no_NA$total_worms_observed
datb_no_NA$fra_adult <- datb_no_NA$Total_adults/datb_no_NA$total_worms_observed

#make variables to put the df's in

ad_coeff_df <- NULL
L4_coeff_df <- NULL


#for each plate, get the logistic model inflection point (median hour) for the timing of the developmental milestone (L4 or adult onset)
for (i in levels(datb_no_NA$unique_plate_id)){
	the_plate <- datb_no_NA[datb_no_NA$unique_plate_id == i,]
	
	#get minimum number of worms observed
	min_worms_plate <- min(the_plate$total_worms_observed)
	
	#normalize to minimum number of worms observed
	
	the_plate$L1_L3_norm <- round(the_plate$fra_L1_L3*min_worms_plate)
	the_plate$L4_older_norm <- round(the_plate$fra_L4_older*min_worms_plate)
	the_plate$total_adult_norm <- round(the_plate$fra_adult*min_worms_plate)
	
	#get the worms not at milestone
	the_plate$not_L4 <- min_worms_plate-the_plate$L4_older_norm
	the_plate$not_adult <- min_worms_plate-the_plate$total_adult_norm
	#na's are 0
	the_plate[is.na(the_plate)] <- 0
	#expand the number of rows according to number of worms at milestone-- each worm gets a row per time observed
	dev_data_expand_at_L4 <- the_plate[rep(1:nrow(the_plate), the_plate$L4_older_norm),]
	dev_data_expand_at_adult <- the_plate[rep(1:nrow(the_plate), the_plate$total_adult_norm),]
	
	#expand the number of rows according to number of worms NOT at milestone-- each worm gets a row per time observed
	
	dev_data_expand_NOT_at_L4 <- the_plate[rep(1:nrow(the_plate), the_plate$not_L4),]
	dev_data_expand_NOT_at_adult <- the_plate[rep(1:nrow(the_plate), the_plate$not_adult),]
	
	
	#add column of ones to hatched worms
	dev_data_expand_at_L4$milestone_status <- rep(1,nrow(dev_data_expand_at_L4))
	dev_data_expand_at_adult$milestone_status <- rep(1,nrow(dev_data_expand_at_adult))
	
	#add column of zeros to not hatched worms
	dev_data_expand_NOT_at_L4$milestone_status <- rep(0,nrow(dev_data_expand_NOT_at_L4))
	dev_data_expand_NOT_at_adult$milestone_status <- rep(0,nrow(dev_data_expand_NOT_at_adult))
	
	
	L4_data_expand <- rbind(dev_data_expand_at_L4, dev_data_expand_NOT_at_L4)
	Adult_data_expand <- rbind(dev_data_expand_at_adult, dev_data_expand_NOT_at_adult)
		
	L4_coeff <- func_log_inf_conf(L4_data_expand,L4_data_expand$Hours_After_Hatch,L4_data_expand$milestone_status)
	
	Adult_coeff <- func_log_inf_conf(Adult_data_expand,Adult_data_expand$Hours_After_Hatch,Adult_data_expand$milestone_status)
	
	L4_df_to_add <- data.frame(Experiment= unique(as.character(L4_data_expand$Experiment)), Observer=unique(as.character(L4_data_expand$Observer)), Species=unique(as.character(L4_data_expand$Species)), Strain=unique(as.character(L4_data_expand$Strain)), unique_plate_id=unique(as.character(L4_data_expand$unique_plate_id)), Plate_type=unique(as.character(L4_data_expand$Plate_type)), Milestone= "L4", midpoint=L4_coeff[1],midpoint_95_CI_low=L4_coeff[2],midpoint_95_CI_high=L4_coeff[3])

	Adult_df_to_add <- data.frame(Experiment= unique(as.character(Adult_data_expand$Experiment)), Observer=unique(as.character(Adult_data_expand$Observer)), Species=unique(as.character(Adult_data_expand$Species)), Strain=unique(as.character(Adult_data_expand$Strain)), unique_plate_id=unique(as.character(Adult_data_expand$unique_plate_id)), Plate_type=unique(as.character(Adult_data_expand$Plate_type)), Milestone= "Adult", midpoint=Adult_coeff[1],midpoint_95_CI_low=Adult_coeff[2],midpoint_95_CI_high=Adult_coeff[3])

	L4_coeff_df <- rbind(L4_coeff_df,L4_df_to_add)

	ad_coeff_df <- rbind(ad_coeff_df,Adult_df_to_add)
}

#woo!

all_estimates <- rbind(L4_coeff_df,ad_coeff_df)

#write.table(all_estimates,"timing_estimates_per_plate.tsv",row.names = FALSE,col.names=TRUE,sep="\t")



#remove experiment A (it was synchronized in a different way-- bleaching v egg lay)

no_exp_a <- all_estimates[all_estimates$Experiment != "A",]

no_exp_a$Experiment <- as.factor(no_exp_a$Experiment)

no_exp_a$Experiment <- droplevels(no_exp_a$Experiment)

#get new group-- strain plus plate type
no_exp_a$Strain.Plate <- paste(no_exp_a$Strain,no_exp_a$Plate_type)
#re subset
L4_ests <- no_exp_a[no_exp_a$Milestone == "L4",]

adult_ests <- no_exp_a[no_exp_a$Milestone == "Adult",]
#cohen's d effect size, all pairwise comparisons
#set seed for reproducibility, bootstrap
set.seed(666)

L4_cohen_d <- as.data.frame(L4_ests %>% cohens_d(midpoint ~ Strain.Plate,ci = TRUE))

L4_cohen_d

#        .y.         group1         group2    effsize n1 n2 conf.low conf.high
#1  midpoint  MQ130 2-4_DHB      MQ130 NGM -2.3643673  7  7    -6.38     -1.18
#2  midpoint  MQ130 2-4_DHB  NKZ35 2-4_DHB -2.4418444  7  7    -4.90     -1.77
#3  midpoint  MQ130 2-4_DHB      NKZ35 NGM -6.3024950  7  7   -16.39     -4.77
#4  midpoint  MQ130 2-4_DHB PD1074 2-4_DHB  1.2979801  7  6     0.31      4.63
#5  midpoint  MQ130 2-4_DHB     PD1074 NGM  0.3182442  7  6    -1.79      1.17
#6  midpoint      MQ130 NGM  NKZ35 2-4_DHB -1.5391527  7  7    -3.61     -0.82
#7  midpoint      MQ130 NGM      NKZ35 NGM -4.8836947  7  7   -11.74     -3.71
#8  midpoint      MQ130 NGM PD1074 2-4_DHB  3.8798658  7  6     2.74      8.83
#9  midpoint      MQ130 NGM     PD1074 NGM  3.0250767  7  6     2.31      7.33
#10 midpoint  NKZ35 2-4_DHB      NKZ35 NGM -1.1509509  7  7    -3.58     -0.11
#11 midpoint  NKZ35 2-4_DHB PD1074 2-4_DHB  2.9786949  7  6     2.36      5.76
#12 midpoint  NKZ35 2-4_DHB     PD1074 NGM  2.6042404  7  6     2.01      4.88
#13 midpoint      NKZ35 NGM PD1074 2-4_DHB  7.3394934  7  6     5.76     16.98
#14 midpoint      NKZ35 NGM     PD1074 NGM  6.8483227  7  6     5.36     17.22
#15 midpoint PD1074 2-4_DHB     PD1074 NGM -1.1262033  6  6    -3.82      0.14
#   magnitude
#1      large
#2      large
#3      large
#4      large
#5      small
#6      large
#7      large
#8      large
#9      large
#10     large
#11     large
#12     large
#13     large
#14     large
#15     large


adult_cohen_d <- as.data.frame(adult_ests %>% cohens_d(midpoint ~ Strain.Plate,ci = TRUE))

adult_cohen_d

#        .y.         group1         group2    effsize n1 n2 conf.low conf.high
#1  midpoint  MQ130 2-4_DHB      MQ130 NGM -3.3108861  7  7    -5.95     -2.37
#2  midpoint  MQ130 2-4_DHB  NKZ35 2-4_DHB -3.5837758  7  7    -7.98     -2.81
#3  midpoint  MQ130 2-4_DHB      NKZ35 NGM -5.9597447  7  7   -14.57     -5.25
#4  midpoint  MQ130 2-4_DHB PD1074 2-4_DHB  3.3618034  7  6     2.12      7.99
#5  midpoint  MQ130 2-4_DHB     PD1074 NGM  2.5387068  7  6     1.69      5.25
#6  midpoint      MQ130 NGM  NKZ35 2-4_DHB -2.1182201  7  7    -4.28     -1.51
#7  midpoint      MQ130 NGM      NKZ35 NGM -4.0608554  7  7   -10.23     -3.59
#8  midpoint      MQ130 NGM PD1074 2-4_DHB  7.1664374  7  6     5.20     15.83
#9  midpoint      MQ130 NGM     PD1074 NGM  5.8522554  7  6     4.73     11.59
#10 midpoint  NKZ35 2-4_DHB      NKZ35 NGM -0.7046609  7  7    -2.50      0.41
#11 midpoint  NKZ35 2-4_DHB PD1074 2-4_DHB  5.2372995  7  6     4.15     11.90
#12 midpoint  NKZ35 2-4_DHB     PD1074 NGM  4.8463854  7  6     3.89      9.58
#13 midpoint      NKZ35 NGM PD1074 2-4_DHB  8.2845162  7  6     7.12     37.83
#14 midpoint      NKZ35 NGM     PD1074 NGM  7.5750334  7  6     6.78     16.87
#15 midpoint PD1074 2-4_DHB     PD1074 NGM -0.5442514  6  6    -2.50      0.67
#   magnitude
#1      large
#2      large
#3      large
#4      large
#5      large
#6      large
#7      large
#8      large
#9      large
#10  moderate
#11     large
#12     large
#13     large
#14     large
#15  moderate


L4_wilcox <- as.data.frame(L4_ests %>% wilcox_test(midpoint ~ Strain.Plate,p.adjust.method="BH"))

adult_wilcox <- as.data.frame(adult_ests %>% wilcox_test(midpoint ~ Strain.Plate,p.adjust.method="BH"))

L4_wilcox

#        .y.         group1         group2 n1 n2 statistic        p p.adj
#1  midpoint  MQ130 2-4_DHB      MQ130 NGM  7  7         3 0.004000 0.006
#2  midpoint  MQ130 2-4_DHB  NKZ35 2-4_DHB  7  7         0 0.000583 0.002
#3  midpoint  MQ130 2-4_DHB      NKZ35 NGM  7  7         0 0.000583 0.002
#4  midpoint  MQ130 2-4_DHB PD1074 2-4_DHB  7  6        34 0.073000 0.092
#5  midpoint  MQ130 2-4_DHB     PD1074 NGM  7  6        19 0.836000 0.836
#6  midpoint      MQ130 NGM  NKZ35 2-4_DHB  7  7         6 0.018000 0.024
#7  midpoint      MQ130 NGM      NKZ35 NGM  7  7         0 0.000583 0.002
#8  midpoint      MQ130 NGM PD1074 2-4_DHB  7  6        42 0.001000 0.002
#9  midpoint      MQ130 NGM     PD1074 NGM  7  6        42 0.001000 0.002
#10 midpoint  NKZ35 2-4_DHB      NKZ35 NGM  7  7        12 0.128000 0.148
#11 midpoint  NKZ35 2-4_DHB PD1074 2-4_DHB  7  6        42 0.001000 0.002
#12 midpoint  NKZ35 2-4_DHB     PD1074 NGM  7  6        42 0.001000 0.002
#13 midpoint      NKZ35 NGM PD1074 2-4_DHB  7  6        42 0.001000 0.002
#14 midpoint      NKZ35 NGM     PD1074 NGM  7  6        42 0.001000 0.002
#15 midpoint PD1074 2-4_DHB     PD1074 NGM  6  6         9 0.180000 0.193
#   p.adj.signif
#1            **
#2            **
#3            **
#4            ns
#5            ns
#6             *
#7            **
#8            **
#9            **
#10           ns
#11           **
#12           **
#13           **
#14           **
#15           ns


adult_wilcox


#        .y.         group1         group2 n1 n2 statistic        p p.adj
#1  midpoint  MQ130 2-4_DHB      MQ130 NGM  7  7         0 0.000583 0.001
#2  midpoint  MQ130 2-4_DHB  NKZ35 2-4_DHB  7  7         0 0.000583 0.001
#3  midpoint  MQ130 2-4_DHB      NKZ35 NGM  7  7         0 0.000583 0.001
#4  midpoint  MQ130 2-4_DHB PD1074 2-4_DHB  7  6        42 0.001000 0.001
#5  midpoint  MQ130 2-4_DHB     PD1074 NGM  7  6        42 0.001000 0.001
#6  midpoint      MQ130 NGM  NKZ35 2-4_DHB  7  7         1 0.001000 0.001
#7  midpoint      MQ130 NGM      NKZ35 NGM  7  7         0 0.000583 0.001
#8  midpoint      MQ130 NGM PD1074 2-4_DHB  7  6        42 0.001000 0.001
#9  midpoint      MQ130 NGM     PD1074 NGM  7  6        42 0.001000 0.001
#10 midpoint  NKZ35 2-4_DHB      NKZ35 NGM  7  7        10 0.073000 0.078
#11 midpoint  NKZ35 2-4_DHB PD1074 2-4_DHB  7  6        42 0.001000 0.001
#12 midpoint  NKZ35 2-4_DHB     PD1074 NGM  7  6        42 0.001000 0.001
#13 midpoint      NKZ35 NGM PD1074 2-4_DHB  7  6        42 0.001000 0.001
#14 midpoint      NKZ35 NGM     PD1074 NGM  7  6        42 0.001000 0.001
#15 midpoint PD1074 2-4_DHB     PD1074 NGM  6  6        13 0.485000 0.485
#   p.adj.signif
#1            **
#2            **
#3            **
#4            **
#5            **
#6            **
#7            **
#8            **
#9            **
#10           ns
#11           **
#12           **
#13           **
#14           **
#15           ns



L4_cohen_d[c(10,1,15),]

#        .y.         group1     group2   effsize n1 n2 conf.low conf.high
#10 midpoint  NKZ35 2-4_DHB  NKZ35 NGM -1.150951  7  7    -3.58     -0.11
#1  midpoint  MQ130 2-4_DHB  MQ130 NGM -2.364367  7  7    -6.38     -1.18
#15 midpoint PD1074 2-4_DHB PD1074 NGM -1.126203  6  6    -3.82      0.14
#   magnitude
#10     large
#1      large
#15     large


L4_plot <- L4_cohen_d[c(10,1,15),]

L4_plot <- cbind(L4_plot,Strain=c("C. elegans WT","C. elegans (clk-1)","C. inopinata WT"))

L4_plot
#        .y.         group1     group2   effsize n1 n2 conf.low conf.high
#10 midpoint  NKZ35 2-4_DHB  NKZ35 NGM -1.150951  7  7    -3.58     -0.11
#1  midpoint  MQ130 2-4_DHB  MQ130 NGM -2.364367  7  7    -6.38     -1.18
#15 midpoint PD1074 2-4_DHB PD1074 NGM -1.126203  6  6    -3.82      0.14
#   magnitude             Strain
#10     large      C. elegans WT
#1      large C. elegans (clk-1)
#15     large    C. inopinata WT


adult_plot <- adult_cohen_d[c(10,1,15),]

adult_plot

#        .y.         group1     group2    effsize n1 n2 conf.low conf.high
#10 midpoint  NKZ35 2-4_DHB  NKZ35 NGM -0.7046609  7  7    -2.50      0.41
#1  midpoint  MQ130 2-4_DHB  MQ130 NGM -3.3108861  7  7    -5.95     -2.37
#15 midpoint PD1074 2-4_DHB PD1074 NGM -0.5442514  6  6    -2.50      0.67
#   magnitude
#10  moderate
#1      large
#15  moderate


adult_plot <- cbind(adult_plot,Strain=c("C. elegans WT","C. elegans (clk-1)","C. inopinata WT"))

cd_plot <- rbind(L4_plot,adult_plot)

L4_plot <- L4_cohen_d[c(10,1,15),]

L4_plot <- cbind(L4_plot,Strain=c("C. elegans WT","C. elegans (clk-1)","C. inopinata WT"))

adult_plot <- adult_cohen_d[c(10,1,15),]

adult_plot <- cbind(adult_plot,Strain=c("C. elegans WT","C. elegans (clk-1)","C. inopinata WT"))

L4_plot <- cbind(L4_plot,stage=c("L4","L4","L4"))

adult_plot <- cbind(adult_plot,stage=c("Adult","Adult","Adult"))


cd_plot_2 <- rbind(L4_plot,adult_plot)

cd_plot_2

#         .y.         group1     group2    effsize n1 n2 conf.low conf.high
#10  midpoint  NKZ35 2-4_DHB  NKZ35 NGM -1.1509509  7  7    -3.58     -0.11
#1   midpoint  MQ130 2-4_DHB  MQ130 NGM -2.3643673  7  7    -6.38     -1.18
#15  midpoint PD1074 2-4_DHB PD1074 NGM -1.1262033  6  6    -3.82      0.14
#101 midpoint  NKZ35 2-4_DHB  NKZ35 NGM -0.7046609  7  7    -2.50      0.41
#11  midpoint  MQ130 2-4_DHB  MQ130 NGM -3.3108861  7  7    -5.95     -2.37
#151 midpoint PD1074 2-4_DHB PD1074 NGM -0.5442514  6  6    -2.50      0.67
#    magnitude             Strain stage
#10      large      C. elegans WT    L4
#1       large C. elegans (clk-1)    L4
#15      large    C. inopinata WT    L4
#101  moderate      C. elegans WT Adult
#11      large C. elegans (clk-1) Adult
#151  moderate    C. inopinata WT Adult


cd_plot_2$effsize.neg <- cd_plot_2$effsize*-1
cd_plot_2$conf.low.neg <- cd_plot_2$conf.low*-1
cd_plot_2$conf.high.neg <- cd_plot_2$conf.high*-1
cd_plot_2

#         .y.         group1     group2    effsize n1 n2 conf.low conf.high
#10  midpoint  NKZ35 2-4_DHB  NKZ35 NGM -1.1509509  7  7    -3.58     -0.11
#1   midpoint  MQ130 2-4_DHB  MQ130 NGM -2.3643673  7  7    -6.38     -1.18
#15  midpoint PD1074 2-4_DHB PD1074 NGM -1.1262033  6  6    -3.82      0.14
#101 midpoint  NKZ35 2-4_DHB  NKZ35 NGM -0.7046609  7  7    -2.50      0.41
#11  midpoint  MQ130 2-4_DHB  MQ130 NGM -3.3108861  7  7    -5.95     -2.37
#151 midpoint PD1074 2-4_DHB PD1074 NGM -0.5442514  6  6    -2.50      0.67
#    magnitude             Strain stage effsize.neg conf.low.neg conf.high.neg
#10      large      C. elegans WT    L4   1.1509509         3.58          0.11
#1       large C. elegans (clk-1)    L4   2.3643673         6.38          1.18
#15      large    C. inopinata WT    L4   1.1262033         3.82         -0.14
#101  moderate      C. elegans WT Adult   0.7046609         2.50         -0.41
#11      large C. elegans (clk-1) Adult   3.3108861         5.95          2.37
#151  moderate    C. inopinata WT Adult   0.5442514         2.50         -0.67



cd_plot_2$stage <- as.factor(cd_plot_2$stage)
cd_plot_2$stage <- factor(cd_plot_2$stage, levels=c("L4","Adult"))

cd_plot_2$Strain <- as.factor(cd_plot_2$Strain)
cd_plot_2$Strain <- as.factor(cd_plot_2$Strain)
cd_plot_2$Strain <- factor(cd_plot_2$Strain, levels=c("C. elegans WT", "C. elegans (clk-1)", "C. inopinata WT"))

no_exp_a$Milestone <- factor(no_exp_a$Milestone, levels=c("L4","Adult"))





no_exp_a$Strain <- as.factor(no_exp_a$Strain)
no_exp_a$Milestone <- as.factor(no_exp_a$Milestone)

no_exp_a$Milestone <- factor(no_exp_a$Milestone, levels=c("L4","Adult"))



levels(no_exp_a$Strain)[levels(no_exp_a$Strain)=="MQ130"] <- "C. elegans (clk-1)"
levels(no_exp_a$Strain)[levels(no_exp_a$Strain)=="PD1074"] <- "C. elegans WT"
levels(no_exp_a$Strain)[levels(no_exp_a$Strain)=="NKZ35"] <- "C. inopinata WT"

no_exp_a$Strain <- factor(no_exp_a$Strain, levels=c("C. elegans WT", "C. elegans (clk-1)", "C. inopinata WT"))

#a single plate for figure 1b


plate1 <- datb_no_NA[datb_no_NA$unique_plate_id == "13",]

#get minimum number of worms observed
min_worms_plate <- min(plate1$total_worms_observed)

#normalize to minimum number of worms observed

plate1$L1_L3_norm <- round(plate1$fra_L1_L3*min_worms_plate)
plate1$L4_older_norm <- round(plate1$fra_L4_older*min_worms_plate)
plate1$total_adult_norm <- round(plate1$fra_adult*min_worms_plate)

#get the worms not at milestone
plate1$not_L4 <- min_worms_plate-plate1$L4_older_norm
plate1$not_adult <- min_worms_plate-plate1$total_adult_norm

plate1[is.na(plate1)] <- 0
#expand the number of rows according to number of worms at milestone-- each worm gets a row per time observed
dev_data_expand_at_L4 <- plate1[rep(1:nrow(plate1), plate1$L4_older_norm),]
dev_data_expand_at_adult <- plate1[rep(1:nrow(plate1), plate1$total_adult_norm),]

#expand the number of rows according to number of worms NOT at milestone-- each worm gets a row per time observed

dev_data_expand_NOT_at_L4 <- plate1[rep(1:nrow(plate1), plate1$not_L4),]
dev_data_expand_NOT_at_adult <- plate1[rep(1:nrow(plate1), plate1$not_adult),]


#add column of ones to hatched worms
dev_data_expand_at_L4$milestone_status <- rep(1,nrow(dev_data_expand_at_L4))
dev_data_expand_at_adult$milestone_status <- rep(1,nrow(dev_data_expand_at_adult))

#add column of zeros to not hatched worms
dev_data_expand_NOT_at_L4$milestone_status <- rep(0,nrow(dev_data_expand_NOT_at_L4))
dev_data_expand_NOT_at_adult$milestone_status <- rep(0,nrow(dev_data_expand_NOT_at_adult))


L4_data_expand <- rbind(dev_data_expand_at_L4, dev_data_expand_NOT_at_L4)





#the clk-1 expression data for figure 1a


clk1dat <- read.table("clk-1_fold_change.tsv", sep="\t", header=T)
clk1dat$stage <- factor(clk1dat$stage, levels=c("L3","L4","Adult"))







a <- ggplot(clk1dat, aes(x = stage, y = log2FoldChange)) + geom_col(fill="#6baed6") + geom_hline(yintercept=0,linetype="dashed") + geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = .2) +theme_cowplot() + scale_y_continuous(limits=c(-2,2), breaks=c(-2,-1,0,1,2)) +xlab("Stage") + ylab(expression(Log[2] ~ Fold ~ Change)) + ggtitle("a") + theme(plot.title = element_text(hjust = -0.05))


b <- ggplot(L4_data_expand, aes(x=Hours_After_Hatch, y=milestone_status)) + geom_point(alpha=0.25, size=0.5, position=position_jitter(height=0.03, width=0.2)) + stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE,size=0.75) + geom_hline(yintercept=0.5, linetype="dashed", colour="black") + geom_vline(xintercept=73.72330,linetype="dotted",colour="black") + theme_cowplot() + xlab("Hours since embryo laid") + ylab("Milestone status") + ggtitle("b") + theme(plot.title = element_text(hjust = -0.05))



c <- ggplot(no_exp_a, aes(x=Strain,y=midpoint)) + stat_summary(aes(group=Plate_type),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, colour="black",position = position_dodge(width = 0.9))  + geom_sina(aes(colour=Plate_type),size=1,alpha=0.5,scale="width") + facet_rep_grid(~Milestone) + theme_cowplot() + theme(strip.background =element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +scale_colour_manual(values=c("#3B0000","#FF0000")) +scale_y_continuous(limits=c(0,140),breaks=seq(from=0,to=140,by=20)) + ylab("Hours since embryo laid")  + ggtitle("c") + theme(plot.title = element_text(hjust = -0.05))



d <- ggplot(cd_plot_2,aes(x=Strain, y=effsize.neg,fill=stage)) + geom_hline(yintercept=0,linetype="dashed") + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin=conf.high.neg, ymax=conf.low.neg), colour="black", width=.1, position=position_dodge(0.9)) + theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_manual(values=c("#E69F00","#56B4E9"))  + ylab("Cohen's d effect size") + xlab ("Strain") + scale_y_continuous(limits=c(-1.0,7.0),breaks=seq(from=-1.0,to=7.0,by=1))  + ggtitle("d") + theme(plot.title = element_text(hjust = -0.05))


(a+b)/(c+d)


ggsave("figure_5-15-24.png",width=10, height=9,units="in")

ggsave("figure_5-15-24.pdf",width=10, height=9,units="in")





#relevant stats for the main text


PD1074dat <- no_exp_a[no_exp_a$Strain == "C. elegans WT",]

PD1074datad <- PD1074dat[PD1074dat$Milestone == "Adult",]

aggregate(midpoint~Plate_type,FUN=summary,data=PD1074datad)

#  Plate_type midpoint.Min. midpoint.1st Qu. midpoint.Median midpoint.Mean
#1    2-4_DHB      56.89222         57.72523        59.07603      60.53466
#2        NGM      56.87244         59.27388        62.18305      63.37081
#  midpoint.3rd Qu. midpoint.Max.
#1         60.73694      69.59813
#2         68.52850      70.03851

1-(59.07603/62.18305)

#0.04996571
	#5% faster
MQ130dat <- no_exp_a[no_exp_a$Strain == "C. elegans (clk-1)",]

MQ130datad <- MQ130dat[MQ130dat$Milestone == "Adult",]

aggregate(midpoint~Plate_type,FUN=summary,data=MQ130datad)

#  Plate_type midpoint.Min. midpoint.1st Qu. midpoint.Median midpoint.Mean
#1    2-4_DHB      70.64228         73.38737        77.68458      77.07827
#2        NGM      86.75642         90.27927        90.78422      92.62502
#  midpoint.3rd Qu. midpoint.Max.
#1         79.93076      84.58481
#2         95.67858      98.91880
1-(77.68458/90.78422)
#[1] 0.1442942
	#14% faster

NKZ35dat <- no_exp_a[no_exp_a$Strain == "C. inopinata WT",]

NKZ35datad <- NKZ35dat[NKZ35dat$Milestone == "Adult",]

aggregate(midpoint~Plate_type,FUN=summary,data=NKZ35datad)


#  Plate_type midpoint.Min. midpoint.1st Qu. midpoint.Median midpoint.Mean
#1    2-4_DHB      97.58003        105.63296       110.54379     114.03881
#2        NGM     115.18361        116.33211       117.53945     122.30915
#  midpoint.3rd Qu. midpoint.Max.
#1        122.96384     132.95425
#2        126.72949     137.31781

1-(110.54379/117.53945)
#[1] 0.05951755
	#6% faster


#Number of plates per group

dat_no_expa <- datb_no_NA[datb_no_NA$Experiment != "A",]

dat_no_expa$Strain.Plate_type <- paste(dat_no_expa$Strain, dat_no_expa$Plate_type)

aggregate(unique_plate_id~Strain.Plate_type,FUN=unique,data=dat_no_expa)

#1     MQ130 2-4_DHB 28, 29, 30, 31, 32, 49, 50
#2         MQ130 NGM 13, 14, 15, 16, 17, 47, 48
#3     NKZ35 2-4_DHB 38, 39, 40, 41, 42, 45, 46
#4         NKZ35 NGM 23, 24, 25, 26, 27, 43, 44
#5    PD1074 2-4_DHB     33, 34, 35, 36, 37, 52
#6        PD1074 NGM     18, 19, 20, 21, 22, 51

#1     MQ130 2-4_DHB 7
#2         MQ130 NGM 7
#3     NKZ35 2-4_DHB 7
#4         NKZ35 NGM 7
#5    PD1074 2-4_DHB 6
#6        PD1074 NGM 6


#number of embryos plated

dat_no_expa_day_zero <- dat_no_expa[dat_no_expa$Hours_After_Hatch == 0,]

summary(dat_no_expa_day_zero$Embryos)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  47.00   61.75   79.00   80.95   89.25  152.00

aggregate(Embryos~Strain.Plate_type,FUN=mean,data=dat_no_expa_day_zero)

#  Strain.Plate_type  Embryos
#1     MQ130 2-4_DHB 80.71429
#2         MQ130 NGM 86.00000
#3     NKZ35 2-4_DHB 64.42857
#4         NKZ35 NGM 82.42857
#5    PD1074 2-4_DHB 93.33333
#6        PD1074 NGM 80.50000



aggregate(Embryos~Strain.Plate_type,FUN=summary,data=dat_no_expa_day_zero)

#  Strain.Plate_type Embryos.Min. Embryos.1st Qu. Embryos.Median Embryos.Mean
#1     MQ130 2-4_DHB     54.00000        72.00000       84.00000     80.71429
#2         MQ130 NGM     61.00000        71.50000       76.00000     86.00000
#3     NKZ35 2-4_DHB     49.00000        53.50000       57.00000     64.42857
#4         NKZ35 NGM     47.00000        58.00000       62.00000     82.42857
#5    PD1074 2-4_DHB     80.00000        83.75000       87.50000     93.33333
#6        PD1074 NGM     53.00000        77.50000       84.00000     80.50000
#  Embryos.3rd Qu. Embryos.Max.
#1        92.00000     99.00000
#2        93.00000    136.00000
#3        75.00000     88.00000
#4       100.00000    152.00000
#5        90.50000    131.00000
#6        88.25000     97.00000

