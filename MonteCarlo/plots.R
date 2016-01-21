library(ggplot2)
library(reshape2)
true_spot_price = 10.450576
vol = 0.25
nb_sim_mc = 1000000
setwd("/Users/ericfourrier/Documents/MNF/mnf_e/data")



# Pair plot generators  ---------------------------------------------------

df_uni = read.csv("uniform_congruence.csv", header = TRUE)
df_uni$Generator = 'Congruentiel'
df_halton = read.csv("uniform_halton.csv", header = TRUE)
df_halton$Generator = 'Halton'
df_sqrt = read.csv("uniform_sqrt.csv", header = TRUE)
df_sqrt$Generator = 'SQRT'
df_points_t = rbind.data.frame(df_uni,df_halton,df_sqrt)
# df_points_t$nb_points_criterion = 1:nrow(df_points_t)
# df[df_points_t$nb_points_criterion<100 ]





#Distribution by generators  ---------------------------------------------

df_random_numbers = read.csv("random_number3.csv", header = TRUE)
# ks.test(df_random_numbers$SQRT, "pnorm")
# ks.test(x, "pnorm")
# ks.test(x, "pnorm")
# shapiro.test(df_random_numbers$Normal)
df_random_numbers_m <- melt(df_random_numbers,variable.name = "type_generator")
shapiro.test(df_random_numbers$Normal[0:1000])
ggplot(df_random_numbers_m,aes(x=value,y=..density..))+
    geom_histogram(fill="cornsilk", colour="grey60", size=.2)+
    geom_density(colour="black", size= .4,kernel="gaussian") + 
    facet_wrap(~type_generator,scales = "free") +
    ggtitle("Histogrammes et estimateurs de densité pour nos génerateurs aléatoires")

ggplot(df_points_t, aes(x=Seed1, y=Seed2, colour =Generator)) +
    geom_point() +
    facet_wrap(facets=~Generator, ncol=2,nrow = 2)


# Shapiro and independence test  ------------------------------------------
shapiro.test(df_random_numbers$Normal[1:5000])
shapiro.test(df_random_numbers$Halton[1:5000])
shapiro.test(df_random_numbers$SQRT[1:5000])
chisq.test(df_random_numbers$Normal[1:5000],df_random_numbers$Normal[5001:10000])


# Accelerated convergence  ------------------------------------------------


df_gaussian = read.csv("mc_gaussian.csv")
df_sqrt= read.csv("mc_sqrt.csv")
df_halton = read.csv("mc_halton.csv")

df_total = cbind.data.frame(df_gaussian$Nb_iterations, df_gaussian$Price, df_sqrt$Price, df_halton$Price)
colnames(df_total) <- c("Nb_iterations",'Congruence', 'SQRT', 'Halton')
df_total$CiUpper = true_spot_price + (vol/sqrt(1:nb_sim_mc))
df_total$CiLower = true_spot_price - (vol/sqrt(1:nb_sim_mc))
df_total_m <- melt(df_total,variable.name = "type_generator",id.vars = 'Nb_iterations',value.name = 'Price_mc')

# Monte carlo with different type of random number generators   
ggplot(df_total_m,aes(x=Nb_iterations, y=Price_mc, color = type_generator))+ 
    geom_hline(yintercept=true_spot_price) + 
    geom_line() +
    scale_x_log10() + 
    ylim(c(10,11))
    # geom_text(aes(0,true_spot_price,label = true_spot_price, vjust = -1))


# Standard deviation and convergence of error 

std_mc <- function(n){
    return(sqrt(1/length(df_gaussian$Price[n]) *sum(df_gaussian$Price[n]-true_spot_price)^2))
}
std_mc_vec <- sapply(seq(from = 10, to = length(df_gaussian$Price), by = 1),std_mc)
df_std = data.frame(1:length(std_mc_vec),std_mc_vec)
colnames(df_std) <- c("Nb_iterations", "Std")
ggplot(df_std,aes(x=Nb_iterations, y=Std))+ 
    geom_line() + 
    scale_x_log10() +
    scale_y_log10()


# Asian and loopback ------------------------------------------------------

nsim = 2000000
# asian 
df_a_1 = read.csv("Convergence Asian 4.csv", header = F)
colnames(df_a_1) <- c('Nb_iterations','Price')
df_a_1$Seed = 'Seed1'

df_a_2 = read.csv("Convergence Asian 2.csv", header = F)
colnames(df_a_2) <- c('Nb_iterations','Price')
df_a_2$Seed = 'Seed2'

df_a_3 = read.csv("Convergence Asian 3.csv", header = F)
colnames(df_a_3) <- c('Nb_iterations','Price')
df_a_3$Seed = 'Seed3'

df_a_total <- rbind(df_a_1[seq(10,nsim,20),],df_a_2[seq(10,nsim,20),],df_a_3[seq(10,nsim,10),])
#df_a_total_reduce = df_a_total[seq(10,nsim,10),1:3] 


ggplot(df_a_total,aes(x=Nb_iterations, y=Price, color=Seed))+ 
    geom_line() +
    geom_hline(yintercept=4.86184) + 
    scale_x_log10()

df_lb_1 = read.csv("Convergence Lookback.csv", header = F)
colnames(df_lb_1) <- c('Nb_iterations','Price')
df_lb_1$Seed = 'Seed1'

df_lb_2 = read.csv("Convergence Lookback 2.csv", header = F)
colnames(df_lb_2) <- c('Nb_iterations','Price')
df_lb_2$Seed = 'Seed2'

df_lb_3 = read.csv("Convergence Lookback 3.csv", header = F)
colnames(df_lb_3) <- c('Nb_iterations','Price')
df_lb_3$Seed = 'Seed3'


df_lb_total <- rbind(df_lb_1[seq(10,nsim,10),], df_lb_2[seq(10,nsim,10),], df_lb_3[seq(10,nsim,10),])


ggplot(df_lb_total,aes(x=Nb_iterations, y=Price, color=Seed))+ 
    geom_line() +
    geom_hline(yintercept=16.1037) + 
    scale_x_log10()
