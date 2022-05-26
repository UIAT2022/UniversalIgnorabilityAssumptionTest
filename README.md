# Ignorability Assumption Test via Explorer Variables for Unmeasured Confounding in Causal Inference

The following repository cantains codes for the algorithms or experiments to reproduce the results in the NeurIPS 2022 paper: Test for Ignorability Assumption via Universal Instrumental Variables for Unmeasured Confounding in Causal Inferece

Repository 
===========

* Codes : contains the codes for DWHT, CLRT, UIAT. Note that the codes of DWHT and CLRT are the modified version of the original code of Guo et al. (2014), which are written by R code. and They are with the simulation settings for the compliance class model case. The code for UIAT is written by python code. Note that uiat.py requires package rpy2, so that **you have to confirm that R is installed in your PC.**

* Experiments : contains the code for all the numerical studies and real data analysis. They are written by python code and in the type of notebooks. In numerical studies, all the scenarios are in the baseline setting and we leave font "changeable" to denote the varying coefficients.



Requirements
=============

* python 3.7
* rpy2 3.4.5 (with R-4.1.1)
* scikit-learn 0.23.2
* numpy 1.19.5
* pandas 1.1.3
* xgboost 1.4.2
* statsmodels 0.12.0
* scipy 1.5.2


example
=============
```python
import uiat

UIAT = uiat.UniversalIgnorabilityAssumptionTest(cause=A,effect=Y,explorer=Z,covariates=C,dtype="continuous",verbose=False) 
#explorer means an UIV or a surrogate for UIV. dtype means the type of data. Now the package provides options for continuous case or mixed case.

pvalue = UIAT.test()

print(pvalue) #H0: Ignorability assumption is satisfied. vs. H1: not H0.
```
