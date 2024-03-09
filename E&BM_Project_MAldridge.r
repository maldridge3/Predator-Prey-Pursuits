#Define constants
rad_cap <- 0.1 #radius where predator attempts a capture
rad_esc <- c(10,10,10) #radius where predator w/ strat 1,2,3 ceases pursuit
tmax = 180 #predator's quit time

columns <- c("pred_type","prey_type","mod_coeff","react_t","d_start","Survive") #names for columns of results matrix
results <- matrix(0,1,length(columns)) #empty results matrix
colnames(results) <- columns


#Loop test values
strat_1 <- c(30,0.25,0.01) #low acceleration, high endurance
strat_2 <- c(12.5,0.3,0.5) #medium acceleration, medium endurance
strat_3 <- c(15,1,3) #high acceleration, high velocity low endurance

strat_mat <- cbind(strat_1,strat_2,strat_3) #matrix of predator values

mod_coeff <- seq(0.6,1.4,0.2) #modifies prey strategy about base strategies
react_t <- seq(-0.5,0.5,0.25) #reaction time (> 0 implies prey-initiated; < 0 implies pred-initiated)
d_start <- seq(1,10,1) #starting distances to try


#Model loop
  for (pred_type in 1:3) {
      for (prey_type in 1:3) {
        for(mod in mod_coeff) { 
            for (r in react_t) {
              for (d in d_start) {
                cat("current conditions:",pred_type,prey_type,r,d,"\n")
                  for (t in seq(0,tmax,0.1)) {
                    x_1_prev <- (strat_mat[1,pred_type]^2*(t-0.1) + strat_mat[1,pred_type]^2*strat_mat[2,pred_type]*(exp(-strat_mat[2,pred_type]*(t-0.1))-1)-strat_mat[1,pred_type]/2*strat_mat[3,pred_type]*(t-0.1)^2)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type])))
                    x_1 <- (strat_mat[1,pred_type]^2*t + strat_mat[1,pred_type]^2*strat_mat[2,pred_type]*(exp(-strat_mat[2,pred_type]*t)-1)-strat_mat[1,pred_type]/2*strat_mat[3,pred_type]*t^2)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type])))
                    x_2 <- d + (strat_mat[1,prey_type]^2*(t-r) + strat_mat[1,prey_type]^2*strat_mat[2,prey_type]*(exp(-strat_mat[2,prey_type]*(t-r))-1)-strat_mat[1,prey_type]/2*strat_mat[3,prey_type]*(t-r)^2)/(strat_mat[1,prey_type] - strat_mat[3,prey_type]/strat_mat[2,prey_type]*(1 + log(strat_mat[1,prey_type]*strat_mat[2,prey_type]/strat_mat[3,prey_type])))
                    if(r>= 0 & t<=r & x_2 > 0){
                      x_2 = d #prey function is parabolic about origin and will sometimes give nonzero values when shifted
                    }
                    dist <- abs(x_2 - x_1)
                    if(x_2 <= x_1){
                      survive <- 0 #predator was significantly faster than prey and passed it up in model
                      break
                    }
                    if(dist <= rad_cap){
                      survive <- 0 #prey was captured
                      break
                     }
                    if(dist >= rad_esc[pred_type]){
                        survive <- 1 #prey escaped
                        break                      
                    }
                    if(t > 0.2 & x_1 < x_1_prev  & (strat_mat[1,pred_type]^2*(1 - exp(-strat_mat[2,pred_type]*t))-strat_mat[1,pred_type]*strat_mat[3,pred_type]*t)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type]))) <= 0.5*strat_mat[1,pred_type]){
                      survive <- 1 #predator got too fatigued and gave up
                      break
                    }
                    if(t == tmax){
                      survive <- 1 #predator gave up; prey escaped
                    }
                    } #end t
                if(results[1,1] == 0){
                  results[1,] <- c(pred_type,prey_type,mod,r,d,survive)
                }
                if(!results[1,1] == 0){
                  results <- rbind(results,c(pred_type,prey_type,mod,r,d,survive))
                }
              } #end starting distance
            } #end reaction time
        } #end mod
      }#end prey type
    } #end pred type

df.results <- as.data.frame(results) #save reuslts matrix as a data.frame in case use ggplot2


#######################
####Analyze results####
#######################

#loop to find which prey was most successful, 2nd most successful for each predator type
surv_rates <- matrix(0,3,3) #empty matrix for number survived of each combination (row = prey)
opt_prey <- c() #empty vector for optimal prey types
sub_prey <- c() #empty vector for second best prey

for (pred_type in 1:3) { 
  for (prey_type in 1:3) {
    surv_rates[prey_type,pred_type] <- sum(df.results[which(df.results[,1]==pred_type & df.results[,2]==prey_type),6])/length(which(df.results[,1]==pred_type & df.results[,2]==prey_type)) #sum the number of each prey type that survive for each pred type and divide by number of combos
  } #end prey_type
  opt_prey[pred_type] <- which.max(surv_rates[,pred_type]) 
  sub_prey[pred_type] <- which(surv_rates[,pred_type] == sort(surv_rates[,pred_type])[2])
} #end pred_type

#for some reason, I get a warning message about replacement here for some parameter values.
  #Everything runs fine, so I don't have a clue why.

#loop to form matrices comparing optimal prey survive to mod, react_t, and d_start
for (j in 1:3) {
  survive <- c()
    for (i in 1:length(mod_coeff)) { #form vector of survival for mod opt
      survive[i] <- sum(df.results[which(df.results[,3] == mod_coeff[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]),6])/length(which(df.results[,3] == mod_coeff[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]))
    } #end i
  pred <- rep(j,length(mod_coeff))
  prey <- rep(opt_prey[j],length(mod_coeff))
    if(j == 1){
      mod_opt_relate1 <- cbind(pred,prey,mod_coeff,survive)
    }
    if(j == 2){
      mod_opt_relate2 <- cbind(pred,prey,mod_coeff,survive)
    }
    if(j == 3){
      mod_opt_relate3 <- cbind(pred,prey,mod_coeff,survive)
    }
  survive <- c()
    for (i in 1:length(react_t)) { #form vector of survival for react opt
      survive[i] <- sum(df.results[which(df.results[,4] == react_t[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]),6])/length(which(df.results[,4] == react_t[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]))
    } #end i
  pred <- rep(j,length(react_t))
  prey <- rep(opt_prey[j],length(react_t))
    if(j == 1){
      react_opt_relate1 <- cbind(pred,prey,react_t,survive)
    }
    if(j == 2){
      react_opt_relate2 <- cbind(pred,prey,react_t,survive)
    }
    if(j == 3){
      react_opt_relate3 <- cbind(pred,prey,react_t,survive)
    }
  survive <- c()
  for (i in 1:length(d_start)) { #form vector of survival for start opt
    survive[i] <- sum(df.results[which(df.results[,5] == d_start[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]),6])/length(which(df.results[,5] == d_start[i] & df.results[,1] == j & df.results[,2] == opt_prey[j]))
  } #end i
  pred <- rep(j,length(d_start))
  prey <- rep(opt_prey[j],length(d_start))
  if(j == 1){
    d_opt_relate1 <- cbind(pred,prey,d_start,survive)
  }
  if(j == 2){
    d_opt_relate2 <-  cbind(pred,prey,d_start,survive)
  }
  if(j == 3){
    d_opt_relate3 <-  cbind(pred,prey,d_start,survive)
  }
  survive <- c()
  for (i in 1:length(mod_coeff)) { #form vector of survival for mod sub
    survive[i] <- sum(df.results[which(df.results[,3] == mod_coeff[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]),6])/length(which(df.results[,3] == mod_coeff[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]))
  } #end i
  pred <- rep(j,length(mod_coeff))
  prey <- rep(sub_prey[j],length(mod_coeff))
  if(j == 1){
    mod_sub_relate1 <- cbind(pred,prey,mod_coeff,survive)
  }
  if(j == 2){
    mod_sub_relate2 <- cbind(pred,prey,mod_coeff,survive)
  }
  if(j == 3){
    mod_sub_relate3 <- cbind(pred,prey,mod_coeff,survive)
  }
  survive <- c()
  for (i in 1:length(react_t)) { #form vector of survival for react sub
    survive[i] <- sum(df.results[which(df.results[,4] == react_t[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]),6])/length(which(df.results[,4] == react_t[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]))
  } #end i
  pred <- rep(j,length(react_t))
  prey <- rep(sub_prey[j],length(react_t))
  if(j == 1){
    react_sub_relate1 <- cbind(pred,prey,react_t,survive)
  }
  if(j == 2){
    react_sub_relate2 <- cbind(pred,prey,react_t,survive)
  }
  if(j == 3){
    react_sub_relate3 <- cbind(pred,prey,react_t,survive)
  }
  survive <- c()
  for (i in 1:length(d_start)) { #form vector of survival for start sub
    survive[i] <- sum(df.results[which(df.results[,5] == d_start[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]),6])/length(which(df.results[,5] == d_start[i] & df.results[,1] == j & df.results[,2] == sub_prey[j]))
  } #end i
  pred <- rep(j,length(d_start))
  prey <- rep(sub_prey[j],length(d_start))
  if(j == 1){
    d_sub_relate1 <- cbind(pred,prey,d_start,survive)
  }
  if(j == 2){
    d_sub_relate2 <-  cbind(pred,prey,d_start,survive)
  }
  if(j == 3){
    d_sub_relate3 <-  cbind(pred,prey,d_start,survive)
  }
} #end j

mod_opt_relate <- as.data.frame(rbind(mod_opt_relate1,mod_opt_relate2,mod_opt_relate3)) #combine matrices for each optimal case for mod,react_t,d_start
react_opt_relate <- as.data.frame(rbind(react_opt_relate1,react_opt_relate2,react_opt_relate3))
d_opt_relate <- as.data.frame(rbind(d_opt_relate1,d_opt_relate2,d_opt_relate3))

mod_sub_relate <- as.data.frame(rbind(mod_sub_relate1,mod_sub_relate2,mod_sub_relate3)) #combine matrices for each 2nd best case for mod,react_t,d_start
react_sub_relate <- as.data.frame(rbind(react_sub_relate1,react_sub_relate2,react_sub_relate3))
d_sub_relate <- as.data.frame(rbind(d_sub_relate1,d_sub_relate2,d_sub_relate3))

library(ggplot2)
library(ggthemes)

mod_opt_plot <- ggplot(mod_opt_relate,aes(mod_coeff,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Mod Coefficient", y = "Ratio Survive", title = "Survival vs Mod Coefficient by Predator Type", subtitle = "(Optimal Prey)", color = "Prey type") +
  theme_clean() 

mod_sub_plot #mod coefficient opt plot

mod_sub_plot <- ggplot(mod_sub_relate,aes(mod_coeff,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Mod Coefficient", y = "Ratio Survive", title = "Survival vs Mod Coefficient by Predator Type", subtitle = "(Sub-Optimal Prey)", color = "Prey type") +
  theme_clean() 

mod_sub_plot #mod coefficient sub plot

react_opt_plot <- ggplot(react_opt_relate,aes(react_t,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Reaction time (r < 0 => pred-init.)", y = "Ratio Survive", title = "Survival vs Reaction Time by Predator Type", subtitle = "(Optimal Prey)", color = "Prey type") +
  theme_clean() 

react_opt_plot #react_t opt plot

react_sub_plot <- ggplot(react_sub_relate,aes(react_t,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Reaction time", y = "Ratio Survive", title = "Survival vs Reaction Time by Predator Type", subtitle = "(Sub-Optimal Prey)", color = "Prey type") +
  theme_clean()

react_sub_plot #react_t sub plot

d_opt_plot <- ggplot(d_opt_relate,aes(d_start,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Starting distance", y = "Ratio Survive", title = "Survival vs Starting Distance by Predator Type", subtitle = "(Optimal Prey)", color = "Prey type") +
  theme_clean() +
  theme()

d_opt_plot #d_start opt plot

d_sub_plot <- ggplot(d_sub_relate,aes(d_start,survive,color = factor(prey))) +
  geom_point() +
  facet_grid(rows = vars(pred))+
  labs(x = "Starting distance", y = "Ratio Survive", title = "Survival vs Starting Distance by Predator Type", subtitle = "(Sub-Optimal Prey)", color = "Prey type") +
  theme_clean() +
  theme()

d_sub_plot #d_start sub plot

cor(d_opt_relate[which(d_opt_relate[,1] == 1, d_opt_relate[,2] == 1),3],d_opt_relate[which(d_opt_relate[,1] == 1, d_opt_relate[,2] == 1),4])
cor(d_opt_relate[which(d_opt_relate[,1] == 3, d_opt_relate[,2] == 1),3],d_opt_relate[which(d_opt_relate[,1] == 3, d_opt_relate[,2] == 1),4])



###################################
######Test loop to find bugs#######
###################################

# pred_type = 3
# prey_type = 3
# t=5
# d=3
# r = -0.5
# mod = 1.2
# for(d in 1:10){
#   cat(d)
# for (t in seq(0,tmax,0.01)) {
#   x_1 <- (strat_mat[1,pred_type]^2*t + strat_mat[1,pred_type]^2*strat_mat[2,pred_type]*(exp(-strat_mat[2,pred_type]*t)-1)-strat_mat[1,pred_type]/2*strat_mat[3,pred_type]*t^2)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type])))
#   x_1_prev <- (strat_mat[1,pred_type]^2*(t-0.1) + strat_mat[1,pred_type]^2*strat_mat[2,pred_type]*(exp(-strat_mat[2,pred_type]*(t-0.1))-1)-strat_mat[1,pred_type]/2*strat_mat[3,pred_type]*(t-0.1)^2)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type])))
#   x_2 <- d + (strat_mat[1,prey_type]^2*(t-r) + strat_mat[1,prey_type]^2*strat_mat[2,prey_type]*mod*(exp(-strat_mat[2,prey_type]*mod*(t-r))-1)-strat_mat[1,prey_type]/2*mod*strat_mat[3,prey_type]*(t-r)^2)/(strat_mat[1,prey_type] - strat_mat[3,prey_type]/strat_mat[2,prey_type]*(1 + log(strat_mat[1,prey_type]*strat_mat[2,prey_type]/strat_mat[3,prey_type])))
#   if(r>0 & t<r & !x_2 == 0){
#     x_2 = d
#   }
#   dist <- abs(x_2 - x_1)
#   cat(x_1,x_2,dist,"\n")
#   if(dist >= rad_esc[pred_type]){
#     cat("Oh, hooray", "\n")
#     break
#   }
#   if(t > 0.2 & x_1 < x_1_prev & (strat_mat[1,pred_type]^2*(1 - exp(-strat_mat[2,pred_type]*t))-strat_mat[1,pred_type]*strat_mat[3,pred_type]*t)/(strat_mat[1,pred_type] - strat_mat[3,pred_type]/strat_mat[2,pred_type]*(1 + log(strat_mat[1,pred_type]*strat_mat[2,pred_type]/strat_mat[3,pred_type]))) <= 0.5*strat_mat[1,pred_type]){
#     cat("Oh, huzzah", "\n")
#     break
#   }
#   if(t == tmax){
#     cat("Oh, absolutely", "\n")
#     break
#   }
#   if(dist <= rad_cap){
#     cat("Dang", "\n")
#     break
#   }
#   if(x_1 >= x_2){
#     cat("Snap", "\n")
#     break
#   }
#   if(x_2 <= 0){
#     cat("Oof", "\n")
#     break
#   }
# } #end t
# } #end d
# 
# x_2


