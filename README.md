# compartmentalized_mfa
Compartmentalized Metabolic flux analysis
The most direct approach for quantifying intracellular metabolic flux is isotope tracing. The method is based on feeding cells with isotopically labeled nutrients, measuring the isotopic labeling of intracellular metabolites, and then computationally inferring flux via metabolic flux analysis (MFA) by predicting the fluxes that generate the best fit to the experimental data. However, utilizing this approach with metabolic measurements performed on a whole-cell level typically limits its applicability to inferring whole-cell level metabolic flux – i.e. average flux through all subcellular organelles. Addressing this challenge, we developed a computational method for inferring cytosolic and mitochondrial specific metabolic fluxes based on whole-cell level measurements of metabolite isotopic labeling.
Extending upon standard MFA methodologies, our method further accounts for simulated compartment specific labeling patterns by de-convoluting whole-cell level metabolite isotopic labeling measurements into cytosolic and mitochondrial-specific simulated labeling patterns (Eq. 1). 
(min)┬(v,α) ∑_(i=1)^N▒((X_i^(WC-EXP)-(〖α_i*X〗_i^(CY-COMP) (v)+(1-α_i )*X_i^(MT-COMP) (v)))/S_(X_i^exp ) )^2
s.t.
Sv=0	 	             	Stoichiometric mass balance			
v_lb≤v≤v_ub	       	Lower and upper bound on net fluxes
0≤α≤1	       	Cytosolic relative pool size of metabolite 

where v and α_i denote the simulated fluxes and cytosolic relative pool size of metabolite i respectively. Simulated metabolite mass-isotopomer distributions of cytosolic and mitochondrial metabolites are denoted by X_i^(CY-COMP) and X_i^(MT-COMP) (v) and are uniquely determined by the forward and backward fluxes [ref]. We computed these mass-isotopomer distributions given all forward/backward fluxes via the Elementary Metabolite Unit (EMU approach) [ref]. Measured whole cell level metabolite mass-isotopomer distribution and its standard deviation are denoted by X_i^(WC-EXP) and S_(X_i^exp ). The objective function minimizes the variance-weighted sum of squared residual of the differences between the measured mass-isotopomer distribution of metabolite i measured on a whole-cell level and a convolution of the simulated mitochondrial and cytosolic mass-isotopomer distributions, considering the relative pool size of metabolite α_i in each compartment. The Stoichiometric matrix is denoted by S.  
The non-convex optimization problem was solved using Matlab’s Sequential Quadratic Optimization (SQP), starting from multiple sets of random fluxes to overcome potential local minima. To compute confidence intervals for the ratio of SHMT1 and SHMT2, SQP was iteratively run to compute the maximum log-likelihood estimation while constraining the fluxes ratio to increasing (and then decreasing) values (with a step size equal to 5% of the fluxes ratio predicted in the initial maximum log-likelihood estimation) [ref]. Confidence interval bounds were determined based on the 95% quantile of χ2-distribution with one degree of freedom.

Isotope tracing was performed by feeding exponentially growing cells with [2,3,3-2H]-serine. Whole-cell level mass isotopomer distribution for Serine, Glycine and TTP were measured using LC-MS. To further constraint the model, uptake/secretion rates for Serine, Glycine and Folate were measured. We constructed a compartmentalized one-carbon metabolism model (Figure 1) and inferred the confidence intervals for the ratio of SHMT1 and SHMT2 under low (200nM) and high (2uM) folate concentrations.


