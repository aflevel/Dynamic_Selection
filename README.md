# Case Studies: Dynamic of Adaptation

# Scenario 1: Consequence of K- vs r-strategy on the fate of genetic diversity under selection
	Some organisms favour an abundant reproduction of a few individual genotypes (r startegy like for most plants and marine wildlife) while some other favour having a lot of breeding individuals having a few offsrping they care more about (K strategy like most mammals).
 What would be the best startegy to preserve genetic diversity under increasingly stringent stabilising selection?
*Keywords: k startegy, r stargegy, genetic diversity, stabilising selection

**QUESTIONS**
Q1- Model increasingly stringent stabilising selection:
	Which parameter(s) should be allowed to vary?
	Which range of variation have you decided to test? (tip: you should select both stabilising selection around the initial population mean:0 and further away:20)

Q2- Based on bibliographic references, Which prediction can you make about the level of genetic diversity for each of the alternative startegy?

Q3- Which statistics of the output can be used to assess genetic diversity?
	Which statistical test have you decided to use to test the difference in outcomes?

Q4- Sensitivity analysis:
  Which other one paramter of the model would you investigate to make sure your conclusion is valid?

# Scenario 2: Effect of two major genes versus a polygenic inheritance on the response to selection
Some traits are controlled by a few genes that have a huge effect on the phenotype, called major genes, these are often involved in mendelian genetics as their effect and location in the genome are easily determined. This is the case for genes confering Huntington's disease for example. Alternatively, more subtle variation can be controlled by multiple genes called polygenic inheritance. This is the case for genes predisposing to colorectal cancer for example.
Which of these genetic architecture can respond the quicker to selection
*Keywords: genetic architecture, major gene, polygenic inheritance, response to selection

	
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

**QUESTIONS**
Q1- Model alternative genetic architecture:
	The challenge is to compare the two genetic architectures (2 loci vs multiple loci involved) keeping the initial heritability equal in both cases and significantly different from 0 or 1.
	How would you parameter the standard deviation of the environmental effect as a function of Loci Number? (tip: look at the parameters given above)
Q2- Based on bibliographic references, which prediction can you make about the efficiency of selection for the two types of genetic architecture?
Q3- Which statistics of the output can be used to assess the response to selection and its efficiency? (tip: think about the variation of the summary statistic over the generations)
	Which statistical test have you decided to use to test the difference in outcomes?
Q4- Sensitivity analysis:
	Which other one parameter of the model would you investigate to make sure your conclusion is valid?

# **Scenario 3: Effect of the different selection regimes on the genetic diversity
For similar levels

**QUESTIONS**
Q1- Model alternative genetic architecture:
	The challenge is to compare the three types of selection regime keeping the fraction of the population selected equal across simulation
	How would you parameter the selection differential in each of the 3 scenarios?

Q2- Based on bibliographic references, which prediction can you make about the effect of the different selection regimes?

Q3- Which statistics of the output can be used to assess the level of genetic diversity in each scenario?
	Which statistical test have you decided to use to test the difference in outcomes?

Q4- Sensitivity analysis:
	Which other one parameter of the model would you investigate to make sure your conclusion is valid?


# The Report:
* the assessment will be based on a written component that will include the response to the questions from the scenario that you have to reshape to fit an article format.
* the report should include:
 	- An introduction stating the biological problem you are investigating (based on further research using keywords from the scenario introductory sentence) and describing the research hypothesis based on Q2
 	- A description of the methods you are using, first stating what are the parameters you are interested in (based on Q1) and which methods you will use to test their effect (based on Q3)
	- A result section where you should report your main finding using at least one original figure (not produced by the script) and one table. You can still use the figure from the output but not only!
	- A discussion section where you discuss the significance of the result (eg: does it match your initial prediction?) and provide support for this result using the sensitivity analysis from Q4.
 
