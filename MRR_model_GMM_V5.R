## All Resuls will be stored in the mod.coef.all variable

#### read libraries
library(gmm)
#rm(list=ls())
#getwe()
# tempdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(tempdir)


#### Important for timestamps in POSIXct
options(digits.secs=3)

#Initialise variables
mod.coef.all = array(NA, dim=c(91,8,24,2)) #Dimension have to been known before! 
p.values <- array(NA, dim = c(91,4,24,2))
colnames(mod.coef.all) <- c("theta", "phi", "roh", "alpha","impl.ba", "share.asyinf", "vol.trades", "num.trades")
colnames(p.values) <- c("p.theta", "p.phi", "p.roh", "p.alpha")

# This variables can be changed
lp.h <- 4
years <- c(2015:2016)
hours <- c(1:24)

##### go through every file
#through each year

for(y in years){
  #Through each hour
  for(h in hours){
    print(paste("h =",h))
    
    ##### read and prepare data
    dir <- paste("C:/Users/user/Desktop/Studium/Master/3_Semester/Trading_Room/R_Code/data_year/",y,"_", h,"_","data_active_hour.txt", sep = "")
    data.all <- read.table(dir, sep=" ")
    data.all$delivery_quarter <- NULL
    data.all$is_active <- NULL
    data.all <- data.all[!is.na(data.all$timestamp),]
    data.all$delivery_date <- as.POSIXct(data.all$delivery_date) 
    data.all$timestamp <- as.POSIXct(data.all$timestamp)
    all.days <- unique(data.all$delivery_date)
    
    #we only need transaction data for the model
    data.all <- subset(data.all, data.all$is_transaction == 1)
    
    #go through every day in the file
    for(d in 1:length(all.days)){ #length(all.days))
      #print(paste("d =",d))
      data <- subset(data.all, as.numeric(data.all$delivery_date) == as.numeric(all.days[d]))
      
      ############# create new data.frame (lowest price when multiple exercises)
      tstamp <- as.numeric(as.POSIXct(data$timestamp))
      tstamp_mod <- unique(as.numeric(as.POSIXct(data$timestamp)))
      data_mod <- matrix(NA, nrow = length(tstamp_mod), ncol = 4)
      colnames(data_mod) <- c("timestamp", "is_buyside", "price", "volume")
      data_mod[,1] <- tstamp_mod
      n <- nrow(data_mod)
      for(i in 1:n){
        sub <- subset(data, tstamp == tstamp_mod[i])
        sub <- sub[min(sub$price)==sub$price,]
        data_mod[i,] <- cbind(tstamp_mod[i], sub$is_buyside[1], sub$price[1], sub$volume[1])
      }
      
      
      ########### Create initiation variable
      init_l0 <- -data_mod[,2]
      init_l0 <- replace(init_l0, init_l0 == 0, 1)
      init_l1 <- c(NA, init_l0[-length(init_l0)]) #lag 1
      data_mod <- cbind(data_mod, init_l0, init_l1)
      
      
      ########### Create price lags
      price_l0 <- data_mod[,3]
      price_l1 <- c(NA, price_l0[-length(price_l0)]) #lag 1
      data_mod <- cbind(data_mod, price_l0, price_l1)
      data_mod <- data_mod[-1,] 
      
      ########## define time period which should be modelled
      ltchange <- as.numeric(as.POSIXct("2015-07-16"))
      if(as.numeric(as.POSIXct(data$delivery_date[1]))>ltchange) lt <- 30 else lt <- 45
      ### Define length of end-period
      lp <- lp.h*60*60 
      start <- (data$delivery_hour[1]-1)*60*60+as.numeric(as.POSIXct(data$delivery_date[1]))-lt*60-lp
      end <- (data$delivery_hour[1]-1)*60*60+as.numeric(as.POSIXct(data$delivery_date[1]))-lt*60
      data_mod <- subset(data_mod, data_mod[,1] >= start & data_mod[,1] <= end)
      
      ## Kontrolle
      # plot(data_mod[,3]~data_mod[,1])
      
      
      ################ Fit GMM
      #Create moment conditions (without lambda and constant!)
      data_gmm <- data_mod[,5:8]
      g <- function(tet, x){
        rho <- tet[1]
        phi <- tet[2]
        the <- tet[3]
        alp <- tet[4]
        ut <- x[,3]-x[,4]-(phi+the)*x[,1]+(phi+rho*the)*x[,2]
        m1 <- x[,1]*x[,2]-x[,1]^2*rho
        m2 <- ut-alp
        m3 <- (ut-alp)*x[,1]
        m4 <- (ut-alp)*x[,2]
        f <- cbind(m1, m2, m3, m4)
        return(f)
      }
      
      #estimate the parameters with gmm() if more observation than 30 are available
      if(nrow(data_gmm)>30) res <- gmm(g, data_gmm, c(rho = 0, phi = 0, the = 0, alp = 0)) #type = "iterative", optfct = c("nlminb"), lower = 0, upper = 10
      sum <- summary(res)
      test <- sum$stest$test[1]
      
      res$coefficients
      
      rho <- res$coefficients[[1]]
      phi <- res$coefficients[[2]]
      theta <- res$coefficients[[3]]
      alpha <- res$coefficients[[4]]
      
      #Estimated/implied Bid-Ask-Spread 
      impl.bidask <- 2*(phi+theta)  
      
      #Share of a.i.
      share.ai <- (2*theta)/impl.bidask
      
      #traded volume
      vol.trades <- sum(data_mod[,4]) 
      
      #number of trades
      num.trades <- nrow(data_mod)
      
      #save all data in one matrix
      mod.coef.all[d,,h,(y-(years[1]-1))] <- c(theta, phi, rho, alpha, impl.bidask, share.ai, vol.trades, num.trades)
      p.values[d,,h,(y-(years[1]-1))] <- c(sum$coefficients[1,4],sum$coefficients[2,4],sum$coefficients[3,4],sum$coefficients[4,4])
    }
  }
}
rm(list=c("rho","phi","alpha","theta","sum","test","init_l0","init_l1","price_l0","price_l1","sub","dir","lp","start","end"))


# summary(res)
# #compare rho with real autocorrelation
# cor(data_gmm[, 1], data_gmm[, 2])
# 
# #plot the model fit
# fit_res <- (phi + theta) * data_gmm[, 1] - (phi + theta * rho) * data_gmm[, 2]
# p_diff <- data_gmm[, 3] - data_gmm[, 4]
# layout(1)
# plot(p_diff, main = "Vergleich Model und real Returns", xlab = "Zeitraum der Endperiode (2 Stunden)")
# lines(fit_res, col = 2)
# 
# #normal regression for comparison
# layout(1)
# lm.fit <- lm(p_diff ~ data_gmm[, 1:2] - 1)
# #summary(lm.fit)
# plot(p_diff, main = "Vergleich Model und real Returns", xlab = "Zeitraum der Endperiode (2 Stunden)")
# lines(lm.fit$fitted.values, col = 3)
# lines(fit_res, col=2)


save.image("data_gmm_V5.RData")


########################## plot ######################

#define parameter
rho <- mod.coef.all[,3,,]
phi <- mod.coef.all[,2,,]
theta <- mod.coef.all[,1,,]
alpha <- mod.coef.all[,4,,]
impl.ba <- mod.coef.all[,5,,]
share.asinf <- mod.coef.all[,6,,]
vol.trades <- mod.coef.all[,7,,]
num.trades <- mod.coef.all[,8,,]
min.rho <- min(rho)
max.rho <- max(rho)
min.phi <- min(phi)
max.phi <- max(phi)
min.theta <- min(theta)
max.theta <- max(theta)
min.alpha <- min(alpha)
max.alpha <- max(alpha)
min.impl.ba <- min(impl.ba)
max.impl.ba <- max(impl.ba)
min.share.asinf <- min(share.asinf)
max.share.asinf <- max(share.asinf)

col.h <- rainbow(24)
#Farbübersicht
png(paste("Farbübersicht_stunden.png", sep = ""), height = 800, width = 1200)
layout(1)
plot(1:24, col = col.h, pch = 19, ylim = c(0,30))
legend("topleft", legend = paste("stunde", 1:24), col = col.h, pch = 19, ncol = 4, cex = 3)
dev.off()

#plot of the coefficients without borders
years <- 1:2
for(y in years){
  png(paste("plot_mrr_coefficients_",2014+y,".png", sep = ""), height = 800, width = 1200)
  par(mfrow = c(2,2))
  plot(rho[,1,y], main = paste("Rho - Year", y), ylim = c(min.rho,max.rho), type = "n")
  for(i in 1:24) points(rho[,i,y], col = col.h[i])
  plot(phi[,1,y], main = paste("Phi- Year", y), ylim = c(min.phi,max.phi), type = "n")
  for(i in 1:24) points(phi[,i,y], col = col.h[i])
  plot(theta[,1,y], main = paste("Theta - Year", y), ylim = c(min.theta,max.theta), type = "n")
  for(i in 1:24) points(theta[,i,y], col = col.h[i])
  plot(alpha[,1,y], main = paste("Alpha - Year", y), ylim = c(min.alpha,max.alpha), type = "n")
  for(i in 1:24) points(alpha[,i,y], col = col.h[i])
  dev.off()
  
  png(paste("plot_mrr_bidask_",2014+y,".png", sep = ""), height = 800, width = 1200)
  par(mfrow = c(2,1))
  plot(impl.ba[,1,y], main = paste("implied bid ask, year", y), ylim = c(min.impl.ba,max.impl.ba), type = "n")
  for(i in 1:24) points(impl.ba[,i,y], col = col.h[i])
  plot(share.asinf[,1,y], main = paste("share of asymmetric information, year", y), ylim = c(min.share.asinf,max.share.asinf), type = "n")
  for(i in 1:24) points(share.asinf[,i,y], col = col.h[i])
  dev.off()
}

#plot of the coefficients with boundaries
years <- 1:2
for(y in years){
  png(paste("plot_mrr_coefficients_",2014+y,"._2.png", sep = ""), height = 800, width = 1200)
  par(mfrow = c(2,2))
  plot(rho[,1,y], main = paste("rho, year", y), ylim = c(0,1), type = "n")
  for(i in 1:24) points(rho[,i,y], col = col.h[i])
  plot(phi[,1,y], main = paste("phi, year", y), ylim = c(-1,1), type = "n")
  for(i in 1:24) points(phi[,i,y], col = col.h[i])
  plot(theta[,1,y], main = paste("theta, year", y), ylim = c(-1,1), type = "n")
  for(i in 1:24) points(theta[,i,y], col = col.h[i])
  plot(alpha[,1,y], main = paste("alpha, year", y), ylim = c(-1,1), type = "n")
  for(i in 1:24) points(alpha[,i,y], col = col.h[i])
  dev.off()
  
  png(paste("plot_mrr_bidask_",2014+y,"_2.png", sep = ""), height = 800, width = 1200)
  par(mfrow = c(2,1))
  plot(impl.ba[,1,y], main = paste("implied bid ask, year", y), ylim = c(-1,5), type = "n")
  for(i in 1:24) points(impl.ba[,i,y], col = col.h[i])
  plot(share.asinf[,1,y], main = paste("share of asymmetric information, year", y), ylim = c(-1,1), type = "n")
  for(i in 1:24) points(share.asinf[,i,y], col = col.h[i])
  dev.off()
}

#plot impl.spread against vol
years <- 1:2
hours <- 1:24
png(paste("plot_mrr_spread_vs_vol",y,".png", sep = ""), height = 800, width = 1200)
par(mfrow = c(4,6))
for(y in years){
  for(h in hours){
    plot(impl.ba[,h,y]~vol.trades[,h,y], main = paste("impl.BA vs. vol - h",h,"- y",y))
    abline(lm(impl.ba[,h,y]~vol.trades[,h,y]))
  }
}
dev.off()

#plot impl.spread against number of trades
years <- 1:2
hours <- 1:24
par(mfrow = c(4,6))
for(y in years){
  for(h in hours){
    plot(impl.ba[,h,y]~num.trades[,h,y], main = paste("impl.BA vs. vol - h",h,"- y",y))
    abline(lm(impl.ba[,h,y]~num.trades[,h,y]))
  }
}
#in one plot
years <- 1:2
hours <- 1:24
col.h <- rainbow(24)
png(paste("plot_mrr_spread_vs_numtrades",y,".png", sep = ""), height = 800, width = 1200)
par(mfrow = c(2,1))
for(y in years){
  plot(impl.ba[,h,y]~num.trades[,h,y], type = "n",
       main = paste("Implied Bid-Ask vs. Number of trades - hour",h,"- year",y), 
       ylab = "Implied bid-ask spread",
       xlab = "Number of trades in the last period",
       xlim = c(min(num.trades), max(num.trades)),
       ylim = c(-1, 5))
  #legend("topright", legend = paste("hour", 1:24), col = col.h, pch = 18, ncol = 5, cex = 0.5, lwd = 2)
  for(h in hours){
    points(impl.ba[,h,y]~num.trades[,h,y], col = col.h[h], pch = 18)
    abline(lm(impl.ba[,h,y]~num.trades[,h,y]), col = col.h[h])
  }
}
dev.off()

