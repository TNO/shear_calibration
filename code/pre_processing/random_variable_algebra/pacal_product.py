"""
Generating pdf and cdf of the product of random variables that are inputs in
reliability analyses.
"""
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pacal as p

matplotlib.use("Agg")
plt.interactive(False)
# p.params.general.parallel = False

f_dir = os.path.dirname(__file__)

# discretization
n_discr = int(1e3)

ID = ""
# ID = "test_"

# -----------------------------------------------------------------------------
# Santa's little helpers
# -----------------------------------------------------------------------------


def gumbel_par(x_mean, x_cov):
    gamma = 0.5772156649015328606065120900824024310421
    x_std = x_mean * x_cov
    x_sigma = np.sqrt(6) / np.pi * x_std
    x_mu = x_mean - x_sigma * gamma
    return x_mu, x_sigma


def lognorm_par(x_mean, x_cov):
    x_std = x_mean * x_cov
    x_sigma = np.sqrt(np.log(x_std ** 2 / x_mean ** 2 + 1))
    x_mu = np.log(x_mean ** 2 / np.sqrt(x_std ** 2 + x_mean ** 2))
    return x_mu, x_sigma


# -----------------------------------------------------------------------------
# Traffic load
# -----------------------------------------------------------------------------
# only to get the probability of non-exceedance for the representative value

# Gumbel (load intensity)
X1_mean = 1.00
X1_cov = 0.075

par = gumbel_par(X1_mean, X1_cov)
X1 = p.GumbelDistr(mu=par[0], sigma=par[1])

# Normal (model uncertainty)
X2_mean = 1.0
X2_cov = 0.142

X2 = p.NormalDistr(mu=X2_mean, sigma=X2_cov * X2_mean)

# Product
PX = X1 * X2

# probability of non-exceedance for the representative value (characteristic)
P_repr = 1 - 1/1000
# P_repr = np.exp(-1/1000)

# corresponding representative value (product)
x_repr = PX.get_piecewise_invcdf(use_interpolated=True)(P_repr)

# probability of non-exceedance to get `x_r` from X1
P_repr_X1 = X1.cdf(x_repr)

# -----------------------------------------------------------------------------
# Wind load
# -----------------------------------------------------------------------------

# Gumbel (qref1) - NOT USED IN THE PRODUCT!!
X1_mean = 1.00
X1_cov = 0.27

par = gumbel_par(X1_mean, X1_cov)
X1 = p.GumbelDistr(mu=par[0], sigma=par[1])

# xx = np.array([1, 1.2, 1.8, 2,3])
# X1.pdf(xx)

# Lognormal (ce)
X2_mean = 1.00
X2_cov = 0.15

par = lognorm_par(X2_mean, X2_cov)
NX2 = p.NormalDistr(mu=par[0], sigma=par[1])

# Gumbel (cpe)
X3_mean = 1.00
X3_cov = 0.20

par = gumbel_par(X3_mean, X3_cov)
X3 = p.GumbelDistr(mu=par[0], sigma=par[1])

# Lognormal (cd)
X4_mean = 1.00
X4_cov = 0.15

par = lognorm_par(X4_mean, X4_cov)
NX4 = p.NormalDistr(mu=par[0], sigma=par[1])


# Lognormal (theta_W)
X5_mean = 1.00
X5_cov = 0.10

par = lognorm_par(X5_mean, X5_cov)
NX5 = p.NormalDistr(mu=par[0], sigma=par[1])

# Product
P = X1 * p.exp(NX2) * X3 * p.exp(NX4)

P.plot()
P.summary()


# crazy range
# xmin        = max(P.mean() - 10*P.std(),0)
xmin = P.mean() - 10 * P.std()
xmax = P.mean() + 40 * P.std()


x_grid = np.linspace(xmin, xmax, num=n_discr)
pdf = P.pdf(x_grid)
cdf = P.cdf(x_grid)

# plt.plot(x_grid, pdf)

M = np.column_stack((x_grid, pdf, cdf))
#
np.savetxt(os.path.join(f_dir, "results", ID + "unit_wind_pacal.txt"), M)

# -----------------------------------------------------------------------------
# Snow load
# -----------------------------------------------------------------------------

# Gumbel (sground1)
X1_mean = 1.00
X1_cov = 0.60

par = gumbel_par(X1_mean, X1_cov)
X1 = p.GumbelDistr(mu=par[0], sigma=par[1])

# xx          = np.array([1, 1.2, 1.8, 2,3])
# X1.pdf(xx)

# Normal (ground-to-roof conversion)
X2_mean = 1.00
X2_cov = 0.15

NX2 = p.NormalDistr(mu=X2_mean, sigma=X2_cov * X2_mean)

# Lognormal (theta_S)
X3_mean = 1.00
X3_cov = 0.10

par = lognorm_par(X3_mean, X3_cov)
NX3 = p.NormalDistr(mu=par[0], sigma=par[1])


# Product
P = X1 * NX2

P.plot()
P.summary()

# crazy range
# xmin        = max(P.mean() - 10*P.std(),0)
xmin = P.mean() - 10 * P.std()
xmax = P.mean() + 40 * P.std()

x_grid = np.linspace(xmin, xmax, num=n_discr)
pdf = P.pdf(x_grid)
cdf = P.cdf(x_grid)

plt.plot(x_grid, pdf)

M = np.column_stack((x_grid, pdf, cdf))
#
np.savetxt(os.path.join(f_dir, "results", ID + "unit_snow_pacal.txt"), M)

plt.show()
