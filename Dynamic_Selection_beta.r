#####################################################################################################################################################
# Script for the MSc Class "Case Studies in Bioinformatics" at the University of Melbourne by Alex Fournier-Level                      #
# This scripts needs R/2.9 or above and !!!!!!!!! REQUIRES TO HAVE THE lme4 PACKAGE FROM CRAN INSTALLED!!!!!!!!!!!                     #
# This function simulates the evolution of a single trait in a haploid population under selection                                      #
# Various types of genetic architecture and types of selection can be simulated (directional, stabilising or disruptive)               #
# and the function returns the histogram of trait evolution and the change in heritability. You can also apply a directional climate change #
#####################################################################################################################################################

# The library to be installed to run the model
library(lme4)
library(compiler) #already in the base package
options(warn=-1) #this is to get rid of a warning message which shows up when plotting the results


Selection.Simulator=function (

	# The variables of the model
	GenoNum=10,          # The number of genotypes in the initial population
	GenoRep=5,           # The number of replicate per genotype in the initial population
	LociNum=1,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=0.01,  # The standard deviation of the environmental effect
	GenerationNum=1,     # The number of generation the model is run for
	Fraction=TRUE,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=Inf,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=-Inf,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=TRUE          # Should the function return plots TRUE or FALSE

	) {

# Removing previous graphs
if (Graph==TRUE) while(length(dev.list())) dev.off(dev.list()[1])
# Initialisation of the parameters generation g=1 and no more than one potential mutant MutantNum=1
g=1
MutantNum=1

# First Model the number of Genotypes and number of Replicate per genotypes in the initial population
Gname=rep(paste("Geno",1:GenoNum,sep=""),GenoRep)

# Then model the genetic effect associated with each locus
if (LociFX=="Equal Loci") AllG=c(-sdGeneticFX,sdGeneticFX) # Either as a bimodal effect aka A allele +, T allele -
if (LociFX=="continuous") AllG=c(-rgamma(1000, shape=1, rate = 1),rgamma(1000, shape=1, rate = 1)) # Or as a continuous effect at each locus (allelic series)
G=sample(AllG,1)
GFreq=sample(seq(from=.1,to=.9,by=0.05),1)
GFX=sample(c(0,G),GenoNum,prob=c(GFreq,1-GFreq),replace=T)
GFX=t(matrix(GFX))
colnames(GFX)=unique(Gname)
if (LociNum>1) for (l in 2:LociNum) {
	G=sample(AllG,1)
	GFreq=sample(seq(from=.1,to=.9,by=0.05),1)
	GFX=rbind(GFX,sample(c(0,G),GenoNum,prob=c(GFreq,1-GFreq),replace=T))
}

# Finally model the random environmental variation, different for each individual
EFX=rnorm(GenoNum*GenoRep,0,sdEnvironmentalFX)
names(EFX)=Gname

# Model the Phenotype Y as the sum of the Genetic effects and the Environmental effects
YFX=vector()
for (i in 1:length(Gname)) YFX=c(YFX,sum(GFX[,which(colnames(GFX)==Gname[i])])+EFX[i])

# Drawing the histogram for the initial population
Ymin=min(YFX)
Ymax=max(YFX)
if (Graph==TRUE) {
	if (GenerationNum<6) par(mfrow=c(1,GenerationNum),mar=c(4,3,3,1))
	if (GenerationNum>5&&GenerationNum<20) par(mfrow=c(4,ceiling(GenerationNum/4)),mar=c(4,3,3,1))
	if (GenerationNum<20) hist(YFX,main=paste("Generation",g),xlab="Trait value",ylab="# of Individuals")
}

# Calculating the initial heritability through variance components estimates
heritability=vector()
Data=data.frame(X=names(YFX),Y=YFX)
test=lmer(Y~1+(1|X),data=Data)
Var.Pheno=attr(VarCorr(test),"sc")^2+VarCorr(test)$'X'[1,1]
heritability=VarCorr(test)$'X'[1,1]/Var.Pheno
Ymean=mean(YFX)
Ysd=sd(YFX)

# Transforming the proportion of Phenotype selection into actual phenotypic values
YFXRange=max(YFX)-min(YFX)
if (Sel.min!=-Inf&&Fraction==TRUE) {
	Frac.min=quantile(YFX,Sel.min)
} else { Frac.min=Sel.min }

if (Sel.max!=Inf&&Fraction==TRUE) {
	Frac.max=quantile(YFX,Sel.max)
} else { Frac.max=Sel.max }

# Run the model for the a certain number of generation (GenerationNum)
for (g in 2:GenerationNum) {

if (GenerationNum==1) break

# Modification of the selected fraction following selection
Frac.min=Frac.min+YFXRange*EnvChangeRate
Frac.max=Frac.max+YFXRange*EnvChangeRate

# Subsetting the fraction of phenotypes
	if (Sel.max>Sel.min) {
		YFX=YFX[which(YFX>Frac.min)]
		YFX=YFX[which(YFX<Frac.max)]
	} else {
		YFX1=YFX[which(YFX>Frac.min)]
		YFX2=YFX[which(YFX<Frac.max)]
		YFX=c(YFX1,YFX2)
	}

# Drawing the fraction of phenotype selected at each generation
	if (GenerationNum<20) {
		if (Graph==TRUE) {
			if (Frac.min!=-Inf) abline(v=Frac.min,lty=2,col=g)
			if (Frac.max!=Inf) abline(v=Frac.max,lty=2,col=g)
		}
	}

# If Nobody survives the cut, 
	Gname=names(YFX)
	if (length(Gname)==0) {
		if (Graph==TRUE) print("EXTINCTION :(")
		heritability=c(heritability,0)
		Ymean=c(Ymean,0)
		Ysd=c(Ysd,0)
		break
	}

# Introducing random mutations
	MutationNum=MuRate*GenoNum*GenoRep*LociNum
	Mutate=round(MutationNum)+rbinom(1,1,prob=c(round(MutationNum-round(MutationNum),8),1-round(MutationNum-round(MutationNum),8)))
	if (Mutate>0) {
		while (MutationNum>0) {
			Mutant=sample(unique(Gname)[order(unique(Gname))],1,prob=table(Gname),replace=T)
			MutationFX=sample(AllG,1)
			GFX=cbind(rbind(GFX,rep(0,ncol(GFX))),c(GFX[,Mutant],MutationFX))
			colnames(GFX)[ncol(GFX)]=paste("Mutant",MutantNum,sep="")
			Gname=c(Gname,paste("Mutant",MutantNum,sep=""))
			MutationNum=MutationNum-1
			MutantNum=MutantNum+1
		}
	}

# Resampling the population...
	Gname=sample(unique(Gname)[order(unique(Gname))],GenoNum*GenoRep,prob=table(Gname),replace=T)
#... with new Environmental FX
	EFX=rnorm(GenoNum*GenoRep,0,sdEnvironmentalFX)
	names(EFX)=Gname

	YFX=vector()
#... but with the same Genetic FX
	for (i in 1:length(Gname)) YFX=c(YFX,sum(GFX[,which(colnames(GFX)==Gname[i])])+EFX[i])
	if (GenerationNum<20&&Graph==TRUE) hist(YFX,main=paste("Generation",g),xlim=c(min(Ymin,YFX),max(Ymax,YFX)),xlab="Trait value",ylab="# of Individuals")
# and recomputing the heritability
	Data=data.frame(X=names(YFX),Y=YFX)
	if (length(unique(Data$X))>1) {
		test=lmer(Y~1+(1|X),data=Data)
		Var.Pheno=attr(VarCorr(test),"sc")^2+VarCorr(test)$'X'[1,1]
		h=VarCorr(test)$'X'[1,1]/Var.Pheno
	} else h=0
	heritability=c(heritability,h)
	Ymean=c(Ymean,mean(YFX))
	Ysd=c(Ysd,sd(YFX))
}
#preparing the output
names(heritability)=paste("Generation",1:length(heritability))
out=rbind(heritability,Ymean,Ysd)
GenoFrequency=table(Gname)
GFXFinal=GFX[,colnames(GFX) %in% Gname]
if (!is.matrix(GFXFinal)) {
GFXFinal.name=round(table(GFXFinal)[2]/length(GFXFinal),2)
GFXFinal=max(abs(GFXFinal))
} else {
GFXFinal.name=apply(GFXFinal,1,function(x) {round(table(abs(x))[2]/nrow(GFXFinal),2)})
GFXFinal=apply(GFXFinal,1,function(x) {names(table(abs(x)))[2]})
}
if (length(Gname)>0) {
	names(GFXFinal)=GFXFinal.name
	GFXFinal=na.omit(GFXFinal)
	GFXFinal.name=names(GFXFinal)
	GFXFinal=as.numeric(GFXFinal)
	names(GFXFinal)=GFXFinal.name
	GFXFinal=GFXFinal[order(GFXFinal,decreasing=T)]
} else {
	GenoFrequency="There's no one left...'"
	GFXFinal="... all is lost!"
}
# and a hellish code to get the plots right...
if (Graph==TRUE&&GenerationNum>1) {
	dev.new()
	par(mar=c(5, 4, 4, 4))
	plot(out[1,],type="line",col=2,ylab="Heritability",ylim=c(0,1),xlab="Generation",main="Trait Evolution")
	par(new=TRUE)
	plot(out[2,], axes=FALSE, ylab='', xlab='', bty='c',ylim=c(min(out[2,]-out[3,]),max(out[2,]+out[3,])))
	axis(4)
	mtext("Mean trait value",4,line=3,col=2)
	ypoly=c(out[2,]+out[3,],rev(out[2,]-out[3,]))
	xpoly=c(1:g,g:1)
	polygon(ypoly~xpoly,col="gray",border=NA)
	par(new=TRUE)
	plot(out[2,], axes=FALSE, pch=16,ylab='', xlab='', bty='c',ylim=c(min(out[2,]-out[3,]),max(out[2,]+out[3,])),col=2)
	par(new=TRUE)
	plot(out[1,],type="line",ylab="Heritability",ylim=c(0,1),xlab="Generation",main="Trait Evolution")
	if (LociFX=="continuous") {
}
if (Graph==TRUE&&LociNum>1&&LociFX=="continuous") {
	dev.new()
	barplot(GFXFinal,names.arg=names(GFXFinal),main="Genetic Effect Sizes",ylab="Effect in unit of sd",xlab="Frequency of derived allele")
	}
}

# Finally the objects to be returned
if (LociFX=="Equal Loci") {
	return(list(Trait_Summary=out,Genotype_Frequency=GenoFrequency))
} else {return(list(Trait_Summary=out,Genotype_Frequency=GenoFrequency,Genetic_Effects=GFXFinal))
}
}

# And a compilation of the function to make it slightly faster
Selection.Simulator=cmpfun(Selection.Simulator)


