import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf

globwarm = pd.read_csv("globwarm.csv")
globwarm.head()

lmod = smf.ols(
    formula='nhtemp ~ wusa + jasper + westgreen + chesapeake + tornetrask + urals + mongolia + tasman', data=globwarm).fit()
lmod.summary()

plt.scatter(lmod.resid.iloc[:-1], lmod.resid.iloc[1:])
plt.axhline(0, alpha=0.5)
plt.axvline(0, alpha=0.5)
plt.show()

np.corrcoef(lmod.resid.iloc[:-1], lmod.resid.iloc[1:])

globwarm = globwarm.dropna()
X = sm.add_constant(globwarm.iloc[:, 1:9])
gmod = sm.GLSAR(globwarm.nhtemp, X, rho=1)
res = gmod.iterative_fit(maxiter=6)
gmod.rho

gmod = sm.GLSAR(globwarm.nhtemp, X, rho=1)
for i in range(6):
    results = gmod.fit()
    print("AR coefficients: {0}".format(gmod.rho))
    rho, sigma = sm.regression.yule_walker(results.resid, order=gmod.order)
    gmod = sm.GLSAR(globwarm.nhtemp, X, rho)

oatvar = pd.read_csv("oatvar.csv", index_col=0)
oatvar['variety'] = oatvar['variety'].astype('category')
oatvar['grams'] = oatvar['yield']
oatvar.head()

mmod = smf.mixedlm("grams ~ variety", oatvar, groups=oatvar['block']).fit()
mmod.summary()

ind = sm.cov_struct.Exchangeable()
gmod = smf.gee("grams ~ variety", "block", oatvar, cov_struct=ind).fit()
gmod.summary()

ind.summary()

fpe = pd.read_csv("fpe.csv", index_col=0)
fpe.head()
