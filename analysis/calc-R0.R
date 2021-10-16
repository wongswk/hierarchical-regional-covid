# Calculate R0's (estimates and credible intervals) implied by MCMC samples
# Needs results saved from "run-multiregion.R" and "run-provincewide.R"

D=5;
k1=.2;
k2=1;
q=0.05;
ud=.1;
ur=.02;
e=ud/(ur+ud);

R=function(f, R0b)
{
  B=R0b/(D+1/k2);
  R0=B*( (e^4*(1-e)*(1-f)^2*k1*k2)/((e*(1/D+q)+1)*(e*k1+1)*(e*k2+1)) + (e*f+1-e)^2/(1/D+q) + (e*k1*(e*f+1-e)^2)/(k2*(e*k1+1)*(e*k2+1)) + (e*(e*f+1-e)^2)/((e*k1+1)*(e*k2+1)) + (e*f+1-e)^2/(k2*(e*k1+1)*(e*k2+1)) + (e^2*k1*(e*f^2+1-e))/((e*k1+1)*(e*k2+1))  )
  
  R0
}


# R0 from BC-wide model
mprov <- readRDS("data-generated/fit-endDec-provincewide.rds")
f2p <- mprov$post$f2[,1]
f3p <- mprov$post$f3[,1]
f4p <- mprov$post$f4[,1]
f5p <- mprov$post$f5[,1]
f6p <- mprov$post$f6[,1]
f7p <- mprov$post$f7[,1]
R0p <- mprov$post$R0

Rf2p <- R(f2p, R0p)
Rf3p <- R(f3p, R0p)
Rf4p <- R(f4p, R0p)
Rf5p <- R(f5p, R0p)
Rf6p <- R(f6p, R0p)
Rf7p <- R(f7p, R0p)
str2 <- paste0("(", paste(round(quantile(Rf2p, c(0.025,0.975)),2), collapse=","), ")")
str3 <- paste0("(", paste(round(quantile(Rf3p, c(0.025,0.975)),2), collapse=","), ")")
str4 <- paste0("(", paste(round(quantile(Rf4p, c(0.025,0.975)),2), collapse=","), ")")
str5 <- paste0("(", paste(round(quantile(Rf5p, c(0.025,0.975)),2), collapse=","), ")")
str6 <- paste0("(", paste(round(quantile(Rf6p, c(0.025,0.975)),2), collapse=","), ")")
str7 <- paste0("(", paste(round(quantile(Rf7p, c(0.025,0.975)),2), collapse=","), ")")
paste0(round(mean(Rf2p),2), " & ", round(mean(Rf3p),2), " & ", round(mean(Rf4p),2), " & ", round(mean(Rf5p),2), " & ", round(mean(Rf6p),2), " & ", round(mean(Rf7p),2))
paste0(str2, " & ", str3, " & ", str4, " & ", str5, " & ", str6, " & ", str7)


# R0 for hierarchical regional model
load("data-generated/fit-endDec-multiregion.rda")

Rf2 <- list()
Rf3 <- list()
Rf4 <- list()
Rf5 <- list()
Rf6 <- list()
Rf7 <- list()

for (i in c(2,5,4,1,3)) {
  f2p <- m7$post$f2[,i]
  f3p <- m7$post$f3[,i]
  f4p <- m7$post$f4[,i]
  f5p <- m7$post$f5[,i]
  f6p <- m7$post$f6[,i]
  f7p <- m7$post$f7[,i]
  R0p <- m7$post$R0
  
  Rf2[[i]] <- R(f2p, R0p)
  Rf3[[i]] <- R(f3p, R0p)
  Rf4[[i]] <- R(f4p, R0p)
  Rf5[[i]] <- R(f5p, R0p)
  Rf6[[i]] <- R(f6p, R0p)
  Rf7[[i]] <- R(f7p, R0p)  
  
  str2 <- paste0("(", paste(round(quantile(Rf2[[i]], c(0.025,0.975)),2), collapse=","), ")")
  str3 <- paste0("(", paste(round(quantile(Rf3[[i]], c(0.025,0.975)),2), collapse=","), ")")
  str4 <- paste0("(", paste(round(quantile(Rf4[[i]], c(0.025,0.975)),2), collapse=","), ")")
  str5 <- paste0("(", paste(round(quantile(Rf5[[i]], c(0.025,0.975)),2), collapse=","), ")")
  str6 <- paste0("(", paste(round(quantile(Rf6[[i]], c(0.025,0.975)),2), collapse=","), ")")
  str7 <- paste0("(", paste(round(quantile(Rf7[[i]], c(0.025,0.975)),2), collapse=","), ")")
  show(paste0(round(mean(Rf2[[i]]),2), " & ", round(mean(Rf3[[i]]),2), " & ", round(mean(Rf4[[i]]),2), " & ", round(mean(Rf5[[i]]),2), " & ", round(mean(Rf6[[i]]),2), " & ", round(mean(Rf7[[i]]),2)))
  show(paste0(str2, " & ", str3, " & ", str4, " & ", str5, " & ", str6, " & ", str7))
    
}

# R0 regional-weighted
pop_size <- c(843375,1225195,297570,795116,1889225)
pop_total <- sum(pop_size)
Rf2w <- Rf3w <- Rf4w <- Rf5w <- Rf6w <- Rf7w <- rep(0, length(Rf2[[1]]))
for (i in 1:5) {
  Rf2w <- Rf2w + Rf2[[i]] * pop_size[i]/pop_total
  Rf3w <- Rf3w + Rf3[[i]] * pop_size[i]/pop_total
  Rf4w <- Rf4w + Rf4[[i]] * pop_size[i]/pop_total
  Rf5w <- Rf5w + Rf5[[i]] * pop_size[i]/pop_total
  Rf6w <- Rf6w + Rf6[[i]] * pop_size[i]/pop_total
  Rf7w <- Rf7w + Rf7[[i]] * pop_size[i]/pop_total
}
str2 <- paste0("(", paste(round(quantile(Rf2w, c(0.025,0.975)),2), collapse=","), ")")
str3 <- paste0("(", paste(round(quantile(Rf3w, c(0.025,0.975)),2), collapse=","), ")")
str4 <- paste0("(", paste(round(quantile(Rf4w, c(0.025,0.975)),2), collapse=","), ")")
str5 <- paste0("(", paste(round(quantile(Rf5w, c(0.025,0.975)),2), collapse=","), ")")
str6 <- paste0("(", paste(round(quantile(Rf6w, c(0.025,0.975)),2), collapse=","), ")")
str7 <- paste0("(", paste(round(quantile(Rf7w, c(0.025,0.975)),2), collapse=","), ")")
show(paste0(round(mean(Rf2w),2), " & ", round(mean(Rf3w),2), " & ", round(mean(Rf4w),2), " & ", round(mean(Rf5w),2), " & ", round(mean(Rf6w),2), " & ", round(mean(Rf7w),2)))
show(paste0(str2, " & ", str3, " & ", str4, " & ", str5, " & ", str6, " & ", str7))




