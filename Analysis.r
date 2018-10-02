######################################################################################################################################
#
#Please find the summary of what was (or would have been presented in the Friday the 5th prac)
#if you feel like doing the exercise at the end and send me your code, I'd have a look at it!
#
######################################################################################################################################

#####ML estimate of heritability lognormal when a few loci are involved and turns normal when LociNum ++

Graph=FALSE
par(mfrow=c(1,5))
heri.esti=list()
for (LociNum in c(1,5,10,50,100)) {

print("____________")
print(paste(LociNum,"genes"))

heri.serie=vector()
for (i in 1:100) heri.serie=c(heri.serie,Selection.Simulator()[[1]][1,1])
hist(heri.serie)

loglik.norm <- function(par, y) {sum(dnorm(y, mean=par[1], sd=sqrt(par[2]), log = TRUE))}
MLest.norm <- optim(c(mean = 0, var = 1), fn = loglik.norm,y = heri.serie, control = list(fnscale = -1,reltol = 1e-16))$par
print(loglik.norm(MLest.norm,heri.serie))
heri.esti[[paste(LociNum)]][["norm"]]=MLest.norm

loglik.gamma <- function(par, y) {sum(dgamma(y,shape=par[1],rate=par[2],log=TRUE))}
MLest.gamma <- optim(c(shape = 1, rate = 1), fn = loglik.gamma,y = heri.serie, control = list(fnscale = -1,reltol = 1e-16))$par
print(loglik.gamma(MLest.gamma,heri.serie))
heri.esti[[paste(LociNum)]][["gamma"]]=c(MLest.gamma["shape"]/MLest.gamma["rate"],MLest.gamma["shape"]/MLest.gamma["rate"]^2)
names(heri.esti[[paste(LociNum)]][["gamma"]])=c("mean","var")

loglik.lognorm <- function(par, y) {sum(dnorm(log(1+y), mean=par[1], sd=sqrt(par[2]), log = TRUE))}
MLest.lognorm <- optim(c(mean = 0, var = 1), fn = loglik.lognorm,y = heri.serie, control = list(fnscale = -1,reltol = 1e-16))$par
print(loglik.lognorm(MLest.lognorm,heri.serie))
heri.esti[[paste(LociNum)]][["lognorm"]]=MLest.lognorm

}

######What is the minimal Mu.Rate to cope with climate change ?

GeneNum.list=list()

for (Mu.Rate in seq(0.001,.05,length=6)) {

print("____________")
print(paste("Mutation Rate =",Mu.Rate))

GeneNum.serie=vector()
for (i in 1:100) GeneNum.serie=c(GeneNum.serie,ncol(Selection.Simulator(MuRate=Mu.Rate)[[1]]))
GeneNum.list[[paste(Mu.Rate)]]=length(GeneNum.serie[GeneNum.serie==12])/length(GeneNum.serie)
}

#For the plot
plot(unlist(GeneNum.list)~as.numeric(names(GeneNum.list)),xlog=TRUE,type="line",col=2)

######what is the probability of having an ancestral genotype in the population after 12 generations?

#Tip to test for the presence of a character string in an element, use the following regular expression:
grep("Geno",names(Selection.Simulator()[[2]]))

######For a precise estimate of heritability, is it better to have many genotypes with low replication or many replication of a few genotypes?

######Finally can you draw the evolution of the frequency of each genotype over time?

#Tip:line 163 >GenoFrequency=table(Gname) could be moved to before the closing bracket in line 159 and GenoFrequency should be changed from a simple vector to a list that takes the generation (g) as name/handle

Ymean=matrix(NA,ncol=20)
for (i in 1:100) Ymean=rbind(Ymean,Selection.Simulator(Sel.max=-1,
	                                                   Sel.min=1,
	                                                   Fraction=F,
	                                                   LociNum=10,
	                                                   sdEnvironmentalFX=1,
	                                                   GenerationNum=20,
	                                                   Graph=F)$Trait_Summary["Ymean",])

plot('n',xlim=c(1,20),ylim=c(min(Ymean,na.rm=T),max(Ymean,na.rm=T)))
for (i in 1:nrow(Ymean)) lines(Ymean[i,],col=i)


##########################################################################################################
#                                                 Case Studies:
##########################################################################################################

#Scenario 1: Consequence of K- vs r-strategy on the fate of genetic diversity under selection
# Some organisms favour an abundant reproduction of a few individual genotypes (r startegy like for most plants and marine wildlife) while some other favour having a lot of breeding individuals having a few offsrping they care more about (K strategy like most mammals)
# What would be the best startegy to preserve genetic diversity under increasingly stringent stabilising selection?
# Keywords: k startegy, r stargegy, genetic diversity, stabilising selection

Selection.Simulator(
	GenoNum=10,          # The number of genotypes in the initial population
	GenoRep=50,           # The number of replicate per genotype in the initial population
	LociNum=10,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=1,  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=F,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=.5,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=-.5,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=T          # Should the function return plots TRUE or FALSE
)

#VS
Selection.Simulator(
	GenoNum=50,          # The number of genotypes in the initial population
	GenoRep=10,           # The number of replicate per genotype in the initial population
	LociNum=10,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=1,  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=F,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=.5,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=-.5,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=T          # Should the function return plots TRUE or FALSE
)

#Q1- Model increasingly stringent stabilising selection:
#	Which parameter(s) should be allowed to vary?
#	Which range of variation have you decided to test? (tip: you should select both stabilising selection around the initial population mean:0 and further away:20)
#
#Q2- Based on bibliographic references, Which prediction can you make about the level of genetic diversity for each of the alternative startegy?
#
#Q3- Which statistics of the output can be used to assess genetic diversity?
#	Which statistical test have you decided to use to test the difference in outcomes?
#
#Q4- Sensitivity analysis:
#	Which other one paramter of the model would you investigate to make sure your conclusion is valid?

#Scenario 2: Effect of two major genes versus a polygenic inheritance on the response to selection
# Some traits are controlled by a few genes that have a huge effect on the phenotype, called major genes, these are often involved in mendelian genetics as their effect and location in the genome are easily determined. This is the case for genes confering Huntington's disease for example. Alternatively, more subtle variation can be controlled by multiple genes called polygenic inheritance. This is the case for genes predisposing to colorectal cancer for example.
# Which of these genetic architectures can respond the quicker to selection?
# Keywords: genetic architecture, major gene, polygenic inheritance, response to selection

H2_maj=c()
for (i in 1:100) H2_maj=c(H2_maj,Selection.Simulator(
	GenoNum=50,          # The number of genotypes in the initial population
	GenoRep=10,           # The number of replicate per genotype in the initial population
	LociNum=2,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=sqrt(1),  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=F,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=Inf,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=0,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=F          # Should the function return plots TRUE or FALSE
)[[1]][1,1])

#VS

H2_poly=c()
for (i in 1:100) H2_poly=c(H2_poly,Selection.Simulator(
	GenoNum=50,          # The number of genotypes in the initial population
	GenoRep=10,           # The number of replicate per genotype in the initial population
	LociNum=10,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=sqrt(5),  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=F,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=Inf,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=0,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=F          # Should the function return plots TRUE or FALSE
)[[1]][1,1])

par(mfrow=c(1,2))
hist(H2_maj)
hist(H2_poly)

#QUESTIONS
#
#Q1- Model alternative genetic architecture:
#	The challenge is to compare the two genetic architectures (2 loci vs multiple loci involved) keeping the initial heritability equal in both cases and significantly different from 0 or 1.
#	How would you parameter the standard deviation of the environmental effect as a function of Loci Number? (tip: look at the parameters given above)
#
#Q2- Based on bibliographic references, which prediction can you make about the efficiency of selection for the two types of genetic architecture?
#
#Q3- Which statistics of the output can be used to assess the response to selection and its efficiency? (tip: think about the variation of the summary statistic over the generations)
#	Which statistical test have you decided to use to test the difference in outcomes?
#
#Q4- Sensitivity analysis:
#	Which other one parameter of the model would you investigate to make sure your conclusion is valid?

#Scenario 3: Effect of the different selection regimes on the genetic diversity
#For similar levels of selection in terms of which fraction of the population is selected for, the resulting level of genetic diversity might differ. We have a strong expectation that balancing/disruptive selection preserve more diversity. However, the effect of how disrupted/far appart the two fractions selected are is not clear (eg: selecting from 0-to-25% and 75-100% vs 25-75% the same as selecting 0-45% and 55%-100% vs 5-95%)
#What is the relationship between level of disruption and increase in persiatnce of genetic diversity?
#Keywords: balancing selection, disruptive selection, genetic diversity

# Stabilising
Selection.Simulator(
	GenoNum=50,          # The number of genotypes in the initial population
	GenoRep=10,           # The number of replicate per genotype in the initial population
	LociNum=10,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=1,  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=T,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=???,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=???,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=T          # Should the function return plots TRUE or FALSE
)

# Disruptive
Selection.Simulator(
	GenoNum=50,          # The number of genotypes in the initial population
	GenoRep=10,           # The number of replicate per genotype in the initial population
	LociNum=10,            # The number of loci affecting the trait
	LociFX="Equal Loci",   # The number of allelic states at on locus, options are "Equal Loci" OR "continuous"
	sdGeneticFX=1 ,       # The standard deviation of each individual locus effect
	sdEnvironmentalFX=sqrt(5),  # The standard deviation of the environmental effect
	GenerationNum=12,     # The number of generation the model is run for
	Fraction=T,           # Whether you want to select a fraction of the population or on absolute values for the trait
	Sel.max=???,         # The maximum phenotypic value in the population (if Fraction==TRUE, it is the proportion of initial population
	Sel.min=???,         # The minimum phenotypic value (can be greater than Sel.max to simulate disruptive selection)
	EnvChangeRate=0,     # % of Change in the selected value (as a % of the range of the initial trait value)
	MuRate=0,      # The mutation rate in number of mutation per locus (do not go too high, R is slow when dealing with long loops...)
	Graph=T          # Should the function return plots TRUE or FALSE
)

#QUESTIONS
#
#Q1- Model alternative genetic architecture:
#	The challenge is to compare the three types of selection regime keeping the fraction of the population selected equal across simulation.
#	How would you parameter the selection differential in each of the 3 scenarios so that you can test various level?
#
#Q2- Based on bibliographic references, which prediction can you make about the effect of the different selection regimes and particularly the level of disruption?
#
#Q3- Which statistics of the output can be used to assess the level of genetic diversity in each scenario?
#	Which statistical test have you decided to use to test the difference in outcomes?
#
#Q4- Sensitivity analysis:
#	Which other one parameter of the model would you investigate to make sure your conclusion is valid?


# The Report:
# the assessment will be based on a written component that will include the response to the questions from the scenario that you have to reshape to fit an article format.
# the report should include:
# 	- An introduction stating the biological problem you are investigating (based on further research using keywords from the scenario introductory sentence) and describing the research hypothesis based on Q2
# 	- A description of the methods you are using, first stating what are the parameters you are interested in (based on Q1) and which methods you will use to test their effect (based on Q3)
#	- A result section where you should report your main finding using at least one original figure (not produced by the script) and one table. You can still use the figure from the output but not only!
#	- A discussion section where you discuss the significance of the result (eg: does it match your initial prediction?) and provide support for this result using the sensitivity analysis from Q4.
# 
