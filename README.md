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

**QUESTIONS**
Q1- Model alternative genetic architecture:
	The challenge is to compare the two genetic architectures (2 loci vs multiple loci involved) keeping the initial heritability equal in both cases and significantly different from 0 or 1.
	How would you parameter the standard deviation of the environmental effect as a function of Loci Number? (tip: look at the parameters given above)
Q2- Based on bibliographic references, which prediction can you make about the efficiency of selection for the two types of genetic architecture?
Q3- Which statistics of the output can be used to assess the response to selection and its efficiency? (tip: think about the variation of the summary statistic over the generations)
	Which statistical test have you decided to use to test the difference in outcomes?
Q4- Sensitivity analysis:
	Which other one parameter of the model would you investigate to make sure your conclusion is valid?


# Scenario 3: Effect of the different selection regimes on the genetic diversity
For similar levels of selection in terms of which fraction of the population is selected for, the resulting level of genetic diversity might differ. We have a strong expectation that balancing/disruptive selection preserve more diversity. However, the effect of how disrupted/far appart the two fractions selected are is not clear (eg: selecting from 0-to-25% and 75-100% vs 25-75% the same as selecting 0-45% and 55%-100% vs 5-95%).
What is the relationship between level of disruption and increase in persiatnce of genetic diversity?
*Keywords: balancing selection, disruptive selection, genetic diversity


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
 
