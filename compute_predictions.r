#Giri Gopalan (Joint work with Saku Vrtilek and Luke Bornn)

#PURPOSE: The purpose of this R code is to estimate probabilities that test systems are either a black hole, pulsar, or nonpulsar, in addition to estimates of the standard errors for these probabilities. The combining rules are as in Rubin (1996).
systems <- c('0614p091','serx1','1608m522','aqlx1','1254m69','1916m053','1701m462','gx3p1','gx5m1','1700m37','0535p262','1900m245','0656m072','0115p634','1636m53','LMCX3','1744m28')
blackhole_prob <- rep(0,length(systems))
bh_se <- rep(0,length(systems))
non_pulsar_prob <- rep(0,length(systems))
np_se <- rep(0,length(systems))
pulsar_prob <- rep(0, length(systems))
p_se <- rep(0,length(systems))
for(i in 1:length(systems))
{
    predict_file = paste("./Validation_Output/",systems[i],".RData",sep="");
    load(predict_file);
    #combined estimate is equivalent to the proportion of prediction samples of each type
    blackhole_prob[i] <- length(which(post_pred_samples == 1))/length(post_pred_samples)
    non_pulsar_prob[i] <- length(which(post_pred_samples == 2))/length(post_pred_samples)
    pulsar_prob[i] <- length(which(post_pred_samples == 3))/length(post_pred_samples)
    N <- length(post_pred_samples[,1])
    m <- length(post_pred_samples[1,])
    bh_probs <- apply(post_pred_samples, 2, function(x) length(which(x ==1))/N)
    np_probs <- apply(post_pred_samples, 2, function(x) length(which(x ==2))/N)
    p_probs <- apply(post_pred_samples, 2, function(x) length(which(x ==3))/N)
    #total variance: within imputation variability + between imputation variability (Neglecting correlations! This can be improved!)
    bh_se[i] <- sqrt(sum(bh_probs*(1-bh_probs)/N)/m+((m+1)/(m*(m-1)))*sum((bh_probs-blackhole_prob[i])^2))
    np_se[i] <- sqrt(sum(np_probs*(1-np_probs)/N)/m+((m+1)/(m*(m-1)))*sum((np_probs-non_pulsar_prob[i])^2))
    p_se[i] <-  sqrt(sum(p_probs*(1-p_probs)/N)/m+((m+1)/(m*(m-1)))*sum((p_probs-pulsar_prob[i])^2))
}
predictions <- data.frame(systems,blackhole_prob,bh_se,non_pulsar_prob,np_se,pulsar_prob,p_se)
