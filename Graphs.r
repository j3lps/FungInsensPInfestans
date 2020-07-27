setwd("C:/Users/Joe/source/repos/carolan-papers/x64/Release/")

# Run all the simulations:
# 0: No pathogen
# 1:5: R3:R7, no fungicide
# 6:10: R3:R7, with fungicide, absolute resistance
# 11:15: R3:R7, with 2 fungicides, absolute resistance
# 16:20: R3:R7, with fungicide, partial resistance
# 21:25: R3:R7, with 2 fungicides, partial resistance
for(i in 0:25) system(paste("Kevin.exe",i),show.output.on.console=FALSE)

plotDayProgress<-function(DF,iYear,variable){
	
	# If iYear == NA, then plot the variable over all time
	if(is.na(iYear)){
		
		# Cut out the required variable, DDay and Year
		D.subset = subset(DF,select=c("Year","DDay",variable))
		Time = D.subset$Year + D.subset$DDay / (max(D.subset$DDay)+1)
		
		plot(D.subset[,3] ~ Time,type="l",ylab=variable,xlab="Time (years)")
		
	} else {
		
		D.subset = subset(DF,Year == iYear, select=c("Year","DDay",variable))
		
		plot(D.subset[,3] ~ D.subset$DDay,type="l",ylab=variable, xlab="Time (degree days)")
		
	}
	
}

plotYearProgress<-function(DF,variable){

	obs = DF[[variable]]
	
	plot(obs ~ DF$Year,type="o",ylab=variable,xlab="Time (years)")
	
}

# Plot a variable
plotDayProgress(iYear=1,variable="Severity")

# Plot the severity of the first year of each simulation
# Variables:
# * TotalLeafArea
# * HealthyArea
# * Severity
# * GF0
par(ask=TRUE)
for(i in 0:25){
	D = read.csv(paste("Day",i,".csv",sep=""))
	plotDayProgress(D,NA,"GF0")
}; par(ask=FALSE)

# Loop over each simulation and plot a variable
par(ask=TRUE)
for(i in 0:25){
	Y = read.csv(paste("Year",i,".csv",sep=""))
	plotYearProgress(Y,variable="Severity")
}; par(ask=FALSE)

# Record the year in which severity increases over 1%
effectiveLife = rep(NA,26)
for(i in 0:25){
	# Read in the year
	Y = read.csv(paste("Year",i,".csv",sep=""))
	effectiveLife[i+1] = min(which(Y$Severity > 0.01))
}

# Figure 3. S with no pathogen and SLI with a pathogen but no control
#for(i in 0:1) system(paste("Kevin.exe",i),show.output.on.console=FALSE)
png("Figure3.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Day0.csv")
plot(D$HealthyArea[D$Year==1]~D$DDay[D$Year==1],type="l",xlab="Time (degree days)",ylab="Area index",col="gray",lty=2,ann=FALSE)
mtext(side=1,line=2,"Time (degree days)")
mtext(side=2,line=2,"Area index")
D<-read.csv("Day1.csv")
lines(D$HealthyArea[D$Year==1]~D$DDay[D$Year==1],lty=1)
lines(D$TotalLatent[D$Year==1]~D$DDay[D$Year==1],lty=3)
lines(D$TotalInf[D$Year==1]~D$DDay[D$Year==1],lty=4)
legend("topleft",col=c(1,"gray",1,1),lty=c(1,2,3,4),c("Healthy","Healthy\n(no disease)","Latent","Infectious"),cex=0.8)
dev.off()

# Figure 4. Severity on R3 and R7 no fungicide, and R3 and R7 with fungicide
#for(i in c(1,5,6,10)) system(paste("Kevin.exe",i),show.output.on.console=FALSE)
png("Figure4.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Day1.csv")
plot(100*D$Severity[D$Year==1]~D$DDay[D$Year==1],type="l",xlab="Time (degree days)",ylab="Severity (%)",ylim=c(0,25),ann=FALSE)
mtext(side=1,line=2,"Time (degree days)")
mtext(side=2,line=2,"Severity (%)")
D<-read.csv("Day5.csv")
lines(100*D$Severity[D$Year==1]~D$DDay[D$Year==1])
D<-read.csv("Day6.csv")
lines(100*D$Severity[D$Year==1]~D$DDay[D$Year==1])
D<-read.csv("Day10.csv")
lines(100*D$Severity[D$Year==1]~D$DDay[D$Year==1])
legend("topleft",lty=1,c("R3; no fung", "R7; no fung", "R3; with fung", "R7; with fung"),cex=0.8)
dev.off()

# Figure 5. Resistance frequency over time
png("Figure5.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Year6.csv")
plot(100*D$Gene0~D$Year,type="l",xlim=c(1,15),xlab="Time (Years)",ylab="Resistance frequency (%)",ann=FALSE)
mtext(side=1,line=2,"Time (Years)")
mtext(side=2,line=2,"Resistance frequency (%)")
D<-read.csv("Year7.csv")
lines(100*D$Gene0~D$Year,lty=2)
D<-read.csv("Year8.csv")
lines(100*D$Gene0~D$Year,lty=3)
D<-read.csv("Year9.csv")
lines(100*D$Gene0~D$Year,lty=4)
D<-read.csv("Year10.csv")
lines(100*D$Gene0~D$Year,lty=5)
legend("topleft",lty=1:5,title="Resistance rating:",c("3","4","5","6","7"),cex=0.8)
dev.off()

# Figure 6. Severity over time
png("Figure6.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Year6.csv")
plot(100*D$Severity~D$Year,type="l",xlim=c(1,15),xlab="Time (Years)",ylab="Severity (%)",ann=FALSE)
mtext(side=1,line=2,"Time (Years)")
mtext(side=2,line=2,"Severity (%)")
D<-read.csv("Year7.csv")
lines(100*D$Severity~D$Year,lty=2)
D<-read.csv("Year8.csv")
lines(100*D$Severity~D$Year,lty=3)
D<-read.csv("Year9.csv")
lines(100*D$Severity~D$Year,lty=4)
D<-read.csv("Year10.csv")
lines(100*D$Severity~D$Year,lty=5)
legend("topleft",lty=1:5,title="Resistance rating:",c("3","4","5","6","7"),cex=0.8)
dev.off()

# Figure 7. Resistance rating ~ Effective life (years)
png("Figure7.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
plot(effectiveLife[7:11]~seq(3,7),ylim=c(0,40),xlab="Resistance rating",type="o",ylab="Effective life (years)",ann=FALSE)
mtext(side=1,line=2,"Resistance rating")
mtext(side=2,line=2,"Effective life (years)")
lines(effectiveLife[12:16]~seq(3,7),type="o",lty=2)
lines(effectiveLife[17:21]~seq(3,7),type="o",lty=3)
lines(effectiveLife[22:26]~seq(3,7),type="o",lty=4)
legend("topleft",lty=c(1,2,3,4),c("AR; Solo","AR; Mix", "PR; Solo", "PR; Mix"),cex=0.8)
dev.off()
