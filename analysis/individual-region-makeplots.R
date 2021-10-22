# Make plots based on individual region and hierarchical fits
# Needs results saved from "run-individual-region.R" and "run-multiregion.R"

load("data-generated/fit-endDec-multiregion.rda")
mindiv <- readRDS("data-generated/fit-endDec-individual-regions.rds")
region_names <- c("Island", "Coastal", "Northern", "Interior", "Fraser")

library(dplyr)
library(ggplot2)

.hist_blue <- RColorBrewer::brewer.pal(6, "Blues")[5]
.start <- lubridate::ymd_hms("2020-03-01 00:00:00")
ts_df <- dplyr::tibble(time = m7$time, time_num = seq_along(m7$time))

R0p <- tibble(f=m7$post$R0, para = 1)

library(latex2exp)
for (region in 1:5) {
  # Hierarchical fit
  f2p <- m7$post$f2[,region]
  f3p <- m7$post$f3[,region]
  f4p <- m7$post$f4[,region]
  f5p <- m7$post$f5[,region]
  f6p <- m7$post$f6[,region]
  f7p <- m7$post$f7[,region]  

  # Individual region fit
  f2 <- mindiv[[region]]$post$f2[,1]
  f3 <- mindiv[[region]]$post$f3[,1]
  f4 <- mindiv[[region]]$post$f4[,1]
  f5 <- mindiv[[region]]$post$f5[,1]
  f6 <- mindiv[[region]]$post$f6[,1]
  f7 <- mindiv[[region]]$post$f7[,1]
  R0 <- tibble(f=mindiv[[region]]$post$R0, para = 1)
  
  ggplot(R0, aes(x=f)) + geom_density(alpha = 0.7, fill = .hist_blue, colour = NA, adjust = 1.25) +
    geom_density(data = R0p, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    geom_histogram(aes(y = ..density..),  alpha = 0.3, colour="white", lwd=0.2) +
    coord_cartesian(xlim = c(1.85,3.5), expand = FALSE) +
    xlab("") +   
    scale_x_continuous(breaks = seq(2, 3.4, 0.2)) + 
    ylab("Density") + facet_wrap(~para, labeller = label_bquote(R[0][b])) + ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/individual-R0b-", region_names[region], ".pdf"), width = 4, height = 2.7)
  
  com1 <-tibble(f = c(f2,f3,f4,f5,f6,f7), para = rep(c(2,3,4,5,6,7), each=length(f2)))
  com2 <- tibble(f = c(f2p,f3p,f4p,f5p,f6p,f7p), para = rep(c(2,3,4,5,6,7), each=length(f2p)))

  ggplot(com1, aes(x=f)) + geom_density(alpha = 0.7, fill = .hist_blue, colour = NA, adjust = 1.25) +
    geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    geom_histogram(aes(y = ..density..),  alpha = 0.3, colour="white", lwd=0.2) +
    coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab("Density") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(f[ .(para)])) +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/individual-f-", region_names[region], ".pdf"), width = 7, height = 2.7)

  sampFrac_1p<-m7$post$sampFrac_1[,region]
  sampFrac_2p<-m7$post$sampFrac_2[,region]
  sampFrac_3p<-m7$post$sampFrac_3[,region]
  sampFrac_4p<-m7$post$sampFrac_4[,region]
    
  sampFrac_1<-mindiv[[region]]$post$sampFrac_1[,1]
  sampFrac_2<-mindiv[[region]]$post$sampFrac_2[,1]
  sampFrac_3<-mindiv[[region]]$post$sampFrac_3[,1]
  sampFrac_4<-mindiv[[region]]$post$sampFrac_4[,1]
  com1 <-tibble(f = c(sampFrac_1,sampFrac_2,sampFrac_3,sampFrac_4), R0 = rep(R0$f,4), para = rep(c(1,2,3,4), each=length(sampFrac_1)))
  com2 <- tibble(f = c(sampFrac_1p,sampFrac_2p,sampFrac_3p,sampFrac_4p), para = rep(c(1,2,3,4), each=length(sampFrac_1p)))
  
  ggplot(com1, aes(x=f)) + geom_density(alpha = 1, fill = .hist_blue, colour = NA, adjust = 1.25) + 
    geom_histogram(aes(y = ..density..), fill = .hist_blue, alpha = 0.3, colour="white", lwd=0) + 
    geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab("Density") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(psi[ .(para)]))  +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/individual-psi_r-", region_names[region], ".pdf"), width = 3.4, height = 2.7)
  
  ggplot(com1, aes(x=f, y=R0)) + geom_point(size = 0.1) + 
    #geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    #coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab(TeX("$R_{0b}$")) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(psi[ .(para)]))  +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/individual-R0-vs-psi_r-", region_names[region], ".pdf"), width = 5.7, height = 2.7)
  
}
