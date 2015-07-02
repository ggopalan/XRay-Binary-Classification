#Giri Gopalan (Joint work with Saku Vrtilek and Luke Bornn)

#library(rgl)

# PURPOSE: The training data set consists of data from 24 binary systems. The number of observations from each of these systems
# varies widely, with some systems heavily overrepresented and some systems underrepresented. Additionally we are
# computationally constrained when using a Gaussian process model because of expensive matrix inversion. Hence we sample
# data from each system with probability inversely proportional to the fraction of the total data that the system consists of. 

#Read in the data 
file_list <- list.files("./Training")
fulldata <- {}
for(i in 1:length(file_list)){
dummy <- read.table(paste("./Training/",file_list[i], sep=""))
fulldata <- rbind(fulldata,dummy)
}
#Add a column for star type
compact <- rep(1,dim(fulldata)[1])
hmbh <- which(fulldata[,2] == 'HMBH')
lmbh <- which(fulldata[,2] == 'LMBH')
bh <- c(hmbh,lmbh)
atoll <- which(fulldata[,2] == 'Atoll')
zsource <- which(fulldata[,2] == 'Zsource')
pulsars <- which(fulldata[,2] == 'Pulsar')
nonpulsar <- c(atoll,zsource)
compact[nonpulsar] <- rep(2,length(nonpulsar))
compact[pulsars] <- rep(3,length(pulsars))
fulldata <- cbind(compact,fulldata)
names(fulldata) <- c("Compact_Star","System","Type","Date","CC3","CC1","CC2")
#Remove apparent outliers
fulldata <- fulldata[fulldata[,5]<1.3 & fulldata[,6]<3 & fulldata[,7]<8,]
#Set an upper bound for number of data points to include in a training data set.
data_size <- .10*dim(fulldata)[1]
#Determine vector of probabilities to sample from, which are inversely proprtional to the system size
systems <- levels(fulldata$System)
probs <- sapply(systems, function(x) 1/length(which(fulldata$System == x)))
probs <- probs/sum(probs)
sampling_probs <- sapply(seq(dim(fulldata)[1]), function(x) probs[fulldata[x,2]])
#Sample each data point with the corresponding probability
samples <- sample(dim(fulldata)[1],size=data_size,replace=FALSE,prob=sampling_probs)
training <- fulldata[samples,]
save.image('dat.RData')
