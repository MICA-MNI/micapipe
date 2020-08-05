.. _theory_page:

(Some) Theory
==============================

The following expandable section covers some of the theoretical and mathematical underpinnings of BrainStat. This is obviously neither an exhaustive nor didactially 


Mass-univariate linear modelling 
-------------------

BrainStat implements element-wise univariate and multivariate linear models, similar to one its predecessor `SurfStat <http://www.math.mcgill.ca/keith/surfstat/>`_ 

Linear models describe a continuous response variable as a function of one or more predictors (which can be continuous or factorial). An example of such a model is  

	Y = b0 + b1*x1 + b2*x2 + ... + e 
	
where the b_i are the parameter estimates, x_i are the variables and e represents the error term, which is assumed to be iid. BrainStat has adopted the straightforward formula nomenclature from SurfStat, in which the above model could be specified as 

	Model = 1 + term(x1) + term(x2) + ... 
	
followed by simple model fitting 
	
	slm = BrainStatLinMod....

Within a specified model, one can then interrogate specific contrasts, i.e. effects of variables (or variable combinations) specified in the model. The respective code for this will be. 

	slm = BrainStatT(slm, contrast) 

Where contrast could be something like x1, -x1 from the above model in the case of continuous predictor variables, such as age.  One could also specify the contrast as x1.level1 - x1.level2 should x be a factorial variable. An example could be that x is a variable indicating sex, then the 


Mixed effects models 
-------------------
BrainStat also incorporates element-wise linear mixed effects models, again leveraging functionality from `SurfStat <http://www.math.mcgill.ca/keith/surfstat/>`_ 

Mixed models allow for the incorporation of fixed and random effects, which can be useful when handling hierarchically organized (e.g. longitudinal data) where correlated measures may exist across data points. Mixed effects models also have an advantage of flexibly handing irregular data and missing data points, making these approached more efficient than e.g. repeated measures ANOVA. A mixed effects model decomposes the total variance into fixed and random effects 

	Y = b0 + b1 * x1 + ... + a1 * z1 + ... + e
	
Where the b_i are the parameters associated to the fixed effects specified in the design matrix X, while a_i are the parmeters associated to the random effects specified in the design matrix Z associated for random effects. In BrainStat, implementation would be equivalent to the formula for the simple linear models, but with the addition of random effects i.e. 

	Model = 1 + term(x1) + term(x2) + random(z1) + ... + I 
	
followed by simple model fitting and contrast estimation as above 
	
	slm = BrainStatLinMod....
	slm = BrainStatT(slm, contrast) 


Correction for multiple comparisons  
-------------------
Several ways for multiple comparisons correction are implemented and can be flexibly assessed with one script: 
	slm = BrainStat_multipleComp(slm, type, vargin)
	 
which allows for 'no', 'bonferroni', 'fdr-peak', 'fdr-cluster', 'rft-peak' , 'rft-cluster', 'tfce' as different options. 'No' correction simply places the t-value in a cumulative t-statistics distribution with a given degrees of freedom, and calculates the p-value. 'Bonferroni' multiplifies the p-value by the number of comparisons. 'fdr-peak' implements the Benjamini-Hochberg procedure. 'fdr-cluster' translates this procedure to topological inference. 'rft-peak' and 'rft-cluster' leverage gaussian random field theory. 'tfce' implements threshold free cluster enhancement, which is a permutation based procedure. 


Multivariate associative techniques  
-------------------
BrainStat furthermore implements multivariate associative techniques, namely canonical correlation analysis and partial least squares, which find the optimal mapping between two sets of high dimensional data (e.g. behavioral batteries and brain imaging measures). 

