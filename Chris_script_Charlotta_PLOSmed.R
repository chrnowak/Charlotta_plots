#######################################################################################################################################
###################### Chris, April 2017, PLOS Med revision 1 #########################################################################

library(ggplot2)
library(gridExtra)
library(grid)
library(MendelianRandomization)

############################
## scatter plot function  ##
############################

plot_scatter2 <- function(data, beta_exp, se_exp, beta_out, se_out, p_out, rsid,
                          pleio1, pleio2,
                          title = "PLOT", x_label = "X", y_label ="Y", 
                          p_cutoff = F , 
                          egger_slope = 0, egger_int = 0, ivw_slope = 0){ 
  if(p_cutoff != F){
    data <- data[which(data[,p_out] < p_cutoff), ]
  }
  
  pleio_data <- data[which(data[,pleio1] == 1 & data[,pleio2] == 0),]
  
  df <- data.frame(x = data[,beta_exp], xmin = data[,beta_exp] - (1.96*data[,se_exp]), xmax = data[,beta_exp] + (1.96*data[,se_exp]),
                   y = data[,beta_out], ymin = data[,beta_out] - (1.96*data[,se_out]), ymax = data[,beta_out] + (1.96*data[,se_out]),
                   snp_name = data[,rsid])
  
  df_2 <- data.frame(x = pleio_data[,beta_exp], xmin = pleio_data[,beta_exp] - (1.96*pleio_data[,se_exp]), xmax = pleio_data[,beta_exp] + (1.96*pleio_data[,se_exp]),
                     y = pleio_data[,beta_out], ymin = pleio_data[,beta_out] - (1.96*pleio_data[,se_out]), ymax = pleio_data[,beta_out] + (1.96*pleio_data[,se_out]))
  
  p <- ggplot(data = df,aes(x = x, y = y)) +
    geom_point() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, width = 0)) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, height = 0)) +
    
    # geom_point(data=df_2, aes(x = x, y = y), size = 4, shape = 8, col = "gray20") +
    geom_point(data = df_2, aes(x = x, y = y), col = "blue3") +
    geom_errorbar(data = df_2, aes(ymin = ymin, ymax = ymax, width = 0), col = "blue3") + 
    geom_errorbarh(data = df_2, aes(xmin = xmin, xmax = xmax, height = 0), col = "blue3") +
    
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    geom_abline(slope = egger_slope, intercept = egger_int, col = "gray40", lty= "dotted") +  ## Charlotta, consider changing colour etc. for IVW
    geom_abline(slope = ivw_slope, intercept = 0, col = "gray40", lty = "dashed") +           ## and Egger lines
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
    labs(y = y_label, x = x_label) +
    xlim(0, 0.15) +
    ylim(-0.2, 0.2)
  return(p)
}

############################
## funnel plot function   ##
############################

plot_funnel <- function(data, beta_exp, se_exp, beta_out, se_out, p_out,
                        pleio1, pleio2, 
                        title = "FUNNEL", vertical = F , egger = F, p_cutoff = F){
  if(p_cutoff != F){
    data <- data[which(data[,p_out] < p_cutoff), ]
  }
  
  pleio_data <- data[which(data[,pleio1] == 1 & data[,pleio2] == 0),]
  
  if(vertical == F){
    vertical <- 0
  } 
  if(egger == F){
    egger <- 0
  }
  
  dfunnel <- data.frame(y = data[,beta_exp] / data[,se_out], x = data[,beta_out] / data[,beta_exp],
                        xmin = (data[,beta_out] / data[,beta_exp]) - (1.96*(data[,se_out] / data[,beta_exp])), 
                        xmax = (data[,beta_out] / data[,beta_exp]) + (1.96*(data[,se_out] / data[,beta_exp])))
  
  dfunnel2 <- data.frame(y = pleio_data[,beta_exp] / pleio_data[,se_out], x = pleio_data[,beta_out] / pleio_data[,beta_exp],
                         xmin = (pleio_data[,beta_out] / pleio_data[,beta_exp]) - (1.96*(pleio_data[,se_out] / pleio_data[,beta_exp])), 
                         xmax = (pleio_data[,beta_out] / pleio_data[,beta_exp]) + (1.96*(pleio_data[,se_out] / pleio_data[,beta_exp])))
  
  p2 <- ggplot(data = dfunnel,aes(x = x, y = abs(y))) + 
    geom_point(size = 1) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
    
    geom_point(data = dfunnel2, aes(x = x, y = y), col = "blue3") +
    geom_errorbarh(data = dfunnel2, aes(xmin = xmin, xmax = xmax, height = 0), col = "blue3") +
    # geom_point(data=dfunnel2, aes(x = x, y = y), size = 4, shape = 8, col = "gray20") +
    
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = vertical, col = "gray40", lty = "dashed") +
    geom_vline(xintercept = egger, col = "gray40", lty = "dotted") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "IV strength",x = "IV estimate")  +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
    xlim(-4.7, 4.7) +
    ylim(0, 3)
  return(p2)
}

## upload tove's tables

setwd("/Users/nowakchr/Desktop/Tove_Charlotta_PLOS_Med/")
tab <- read.csv("LF_ex_all_170423_TF.csv", header = T, stringsAsFactors = F)
tab <- tab[-which(tab$t1d_SNP == ""), ]  ## data cleaning - tove's excel file has last column sums

######  Tove: 4 scores: Felix_cons, Felix_liberal, Full_cons, Felix_liberal. generate for both Felix and Full  (in total 4 plots):
######  scatter plot, mark SNPs that are in liberal but not in cons 
######  funnel plot, mark SNPs that are in liberal but not in cons

tab_felix <- tab[which(tab$Felix_liberal == 1),]   # felix dataset
tab_full <- tab[which(tab$Full_liberal == 1),]     # full dataset

## Charlotta's result lists the MR results for Felix, not for FULL hence, for plot, calculate Egger and IVW estimate

mr1 <- mr_allmethods(mr_input(bx = tab_full$ob_Beta, bxse = tab_full$ob_SE, 
                              by = tab_full$t1d_beta, byse = tab_full$t1d_SE))$Values

## note: checking Charlotta's results with FELIX produces the SAME results as in her document :D
## mr_allmethods(mr_input(bx = tab_felix$ob_Beta, bxse = tab_felix$ob_SE, by = tab_felix$t1d_beta, byse = tab_felix$t1d_SE))$Values

## PLOT function below should be clear, 
## p_cutoff         can be set to plot only out-associated SNP > P
## pleio1, pleio2   to mark SNPs in pleio1 but NOT IN pleio2
## remove title etc. put ""

p1 <- plot_scatter2(tab_felix, beta_exp = "ob_Beta", se_exp ="ob_SE", 
                    beta_out = "t1d_beta" , se_out = "t1d_SE", p_out = "t1d_P", rsid = "t1d_SNP", 
                    pleio1 = "Felix_liberal",  pleio2 = "Felix_cons",
                    title ="Felix" , x_label = "Genetic association with BMI", 
                    y_label = "Genetic association with T1D",
                    egger_int = -0.0182446  , egger_slope = 0.7234641  , ivw_slope =   0.4385517) # from Charlotta's results

p2 <- plot_scatter2(tab_full, beta_exp = "ob_Beta", se_exp ="ob_SE",
                    beta_out = "t1d_beta", se_out = "t1d_SE", p_out = "t1d_P", rsid = "t1d_SNP", 
                    pleio1 = "Full_liberal", pleio2 = "Full_cons",
                    title ="Full" , x_label = "Genetic association with BMI", 
                    y_label = "Genetic association with T1D",
                    egger_int = mr1[9,2], egger_slope = mr1[8,2], ivw_slope = mr1[4,2])

p3 <- plot_funnel(tab_felix, beta_exp = "ob_Beta", se_exp ="ob_SE", beta_out = "t1d_beta" ,se_out = "t1d_SE", p_out = "t1d_P",
                  pleio1 = "Felix_liberal", pleio2 = "Felix_cons", 
                  title="",
                  egger = 0.7234641, vertical = 0.4385517)

p4<-plot_funnel(tab_full, beta_exp = "ob_Beta", se_exp ="ob_SE", beta_out = "t1d_beta" , se_out = "t1d_SE", p_out = "t1d_P",
                pleio1 = "Full_liberal", pleio2 = "Full_cons",
                title="",
                egger = mr1[8,2], vertical = mr1[4,2])

# pdf("plot_1.pdf", w = 9, h = 9)
grid.arrange(p1, p2, p3, p4, ncol = 2) #  9 X 9
# dev.off()

## I usually export as PDF to keep files small - for submission, I produce TIFF 300 dpi files by
## opening PDF in Preview (Mac) > Export... > TIFF

## as EPS (looks a bit different)
# postscript("FileName.eps", width = 9, height = 9)
# grid.arrange(p1,p2,p3,p4, ncol=2) #  9 X 9
# dev.off()
            

## make 1 forest plot the Full dataset for  (in total 3 plots)
## 1. Birthweight, mark those SNPs that are in liberal but not in cons
## 2. Education, mark those SNPs that are in liberal but not in cons
## 3. Smoking (OR scale), mark those SNPs that are in liberal but not in cons

## 1 summarizing forest plot such as in the submitted version, 
## but including (in total 1 plot) Felix_cons, Felix_liberal, Full_cons, Full_liberal

############################
## forest plot function   ##
############################

FOREST <- function(data, beta_exp, se_exp, rsid = "t1d_SNP",
                   pleio1, pleio2,
                   title = "PLOT", x_label = "X", y_label ="Y",
                   log = F, hj = 0){ 
  
  data <- data[order(data[, rsid]), ] # alph-num order of rsids
  
  pleio_data <- data[which(data[, pleio1] == 1 & data[, pleio2] == 0),]
  df <- data.frame(x = data[,rsid], y = data[,beta_exp], 
                   ymin = data[,beta_exp] - (1.96*data[,se_exp]), 
                   ymax = data[,beta_exp] + (1.96*data[,se_exp]))
  df_2 <- data.frame(x = pleio_data[,rsid], y = pleio_data[,beta_exp], 
                     ymin = pleio_data[,beta_exp] - (1.96*pleio_data[,se_exp]), 
                     ymax = pleio_data[,beta_exp] + (1.96*pleio_data[,se_exp]))
  x1 <- -0.04
  x2 <- 0.04
  x3 <- 0.033
  i <- 0
  if(log==T){
    df <- data.frame(x = data[,rsid], y = exp(data[,beta_exp]), 
                     ymin = exp(data[,beta_exp] - (1.96*data[,se_exp])), 
                     ymax = exp(data[,beta_exp] + (1.96*data[,se_exp])))
    
    df_2 <- data.frame(x = pleio_data[,rsid],  y = exp(pleio_data[,beta_exp]), 
                       ymin = exp(pleio_data[,beta_exp] - (1.96*pleio_data[,se_exp])), 
                       ymax = exp(pleio_data[,beta_exp] + (1.96*pleio_data[,se_exp])))
    x1 <- 0.87
    x2 <- 1.13
    x3 <- 1.11
    i <- 1
  }
  
  p <- ggplot(data = df, aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
    geom_pointrange(data = df, shape = 18, col ="black") + 
    geom_pointrange(data = df_2, aes(x = x, y = y, ymin = ymin, ymax = ymax), col = "blue3", shape = 18) +
    geom_hline(yintercept=i, lty = 2) +  
    coord_flip() +  
    xlab(y_label) + ylab(x_label) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(x1, x2) +
    geom_text(aes(y = x3, label = 
                    paste(format(round(df$y, digits = 3), nsmall = 3), " [",
                          format(round(df$ymin, digits = 3), nsmall = 3), ", ",
                          format(round(df$ymax, digits = 3), nsmall = 3), "]", sep = "")),
              hjust = hj, size = 3)
  return(p)
}

## PDF 8 x 5
FOREST(tab_full, beta_exp = "bw_beta", se_exp = "bw_SE",
       pleio1 = "Full_liberal", pleio2 = "Full_cons",
       title = "Birthweight",
       x_label = "Beta [95% CI]",
       y_label = "T1D SNP", hj =0.25, log = F)

FOREST(tab_full, beta_exp = "edu_beta", se_exp = "edu_SE",
       pleio1 = "Full_liberal", pleio2 = "Full_cons",
       title = "Educational attainment",
       x_label = "Beta [95% CI]",
       y_label = "T1D SNP", hj =0.3, log = F)

FOREST(tab_full, beta_exp = "smok_beta", se_exp = "smok_SE",
       pleio1 = "Full_liberal", pleio2 = "Full_cons",
       title = "Smoking",
       x_label = "Beta [95% CI]",
       y_label = "T1D SNP", hj =0.3, log = T)


######  Reviewers:
## a) 300 dpi TIF / EPS
## b) Fig 3 and 4:  indicate the pleiotropic variant with different plotting symbol not arrow. 
##                  axis labels could be more descriptive.
## c) Figure 5 - the x-axis should be logarithmic
## d) The figures are somewhat  uninformative (and boring).

######  Charlotta:
##  columns with "tove" ahead is raw data and NOT aligned. T1D_EAF is flipped correctly
##  SNiPA_MAF is MAF for SNIPA and NOT aligned; detailed information about the columns is in Results_170423.
##  F_extended_t1d_170423: 7 Felix SNPs with P  5x10^-8 to 5x10^-6; originally 9 but two removed becaus Locke-P >0.05/9: rs25832 and rs7869969

