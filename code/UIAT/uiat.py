
"""
Universal Ignorability Assumption Test

"""

import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri
from rpy2.robjects.vectors import StrVector
from statsmodels.stats.anova import anova_lm 
import statsmodels.formula.api as smf

required_packages = ['base', 'forecast', 'CondIndTests','devtools'] # list of required R packages 
if all(rpackages.isinstalled(x) for x in required_packages):
    check_packages = True # True if packages are already installed 
else:
   check_packages = False # False if packages are not installed 
if check_packages == False: # Not installed? Then install.
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    packages_to_install = [x for x in required_packages if not rpackages.isinstalled(x)]
    if len(packages_to_install) > 0:
        utils.install_packages(StrVector(packages_to_install))
    check_packages = True 

r = ro.r
base = importr('base')
forecast = importr('forecast')
graphics = importr('graphics')
grdevices = importr('grDevices')
condindtests = importr('CondIndTests')
devtools = importr('devtools')

required_packages = ['RCIT'] # list of required R packages 
if all(rpackages.isinstalled(x) for x in required_packages):
    check_packages = True # True if packages are already installed 
else:
   check_packages = False # False if packages are not installed 
if check_packages == False: # Not installed? Then install.
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    packages_to_install = [x for x in required_packages if not rpackages.isinstalled(x)]
    if len(packages_to_install) > 0:
        r.install_github("ericstrobl/RCIT")
    check_packages = True 

RCIT = importr('RCIT')

class UniversalIgnorabilityAssumptionTest:        #inputs (cause, effect, explorer, covariates) are arrays, dtypes are "continuous" or "mixed".
    def __init__(self, cause=None, effect=None, explorer=None, covariates=None, level=0.05, dtype=None, verbose=True):
        self.A = pd.DataFrame({'A': cause})
        self.Y = pd.DataFrame({'Y': effect})
        self.Z = pd.DataFrame({'Z': explorer})
        self.X = pd.DataFrame(covariates)
        self.level = 0.05
        self.dtype = dtype
        self.verbose = verbose

        
    def test(self): 

        dim_covariates = self.X.shape[1]
        
        if self.dtype=="continuous":
            if dim_covariates==0:
                tmp_conditioning = self.A
            else:
                tmp_conditioning=pd.concat([self.A,self.X],axis=1)
    
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_from_pd_conditioning = ro.conversion.py2rpy(tmp_conditioning)
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_from_pd_effect = ro.conversion.py2rpy(self.Y)
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_from_pd_explorer = ro.conversion.py2rpy(self.Z)
                
            sample_size = r.unlist(r.dim(r_from_pd_conditioning))[0]
            dim = r.unlist(r.dim(r_from_pd_conditioning))[1]
                
            r_from_pd_explorer=r.unlist(r_from_pd_explorer)    
            r_from_pd_effect=r.unlist(r_from_pd_effect)   
            r_from_pd_conditioning=r.unlist(r_from_pd_conditioning)   
            r_from_pd_conditioning=r.matrix(r_from_pd_conditioning,sample_size,dim)
            
            
            test=r.RCoT(r_from_pd_explorer,r_from_pd_effect,r_from_pd_conditioning,seed=2)
            pvalue=r.unlist(test)[0]
            
            if self.verbose == True:
                print(pvalue)
                
        elif self.dtype=="mixed":
            dat_mixed=pd.concat([self.Y,self.Z,self.A,self.X],axis=1)
            dat_mixed.columns = ['col1','col2','col3','col4']
            
            mod_restriced = smf.ols('col1 ~ col3*col4 ', data = dat_mixed).fit()
            mod_unrestriced = smf.ols('col1 ~ col2*col3*col4', data = dat_mixed).fit()
            anovaResults = anova_lm(mod_restriced, mod_unrestriced)
            pvalue=anovaResults.iloc[1,5]
            
            if self.verbose == True:
                print(pvalue)
    
        else:
            raise ValueError("dtype is not valid.")
        
        if self.verbose == True:
            if(pvalue<=self.level):
                print("The Ignorability Assumption is not satisfied.")
            else:
                print("The Ignorability Assumption is satisfied.")
    
        return pvalue
    