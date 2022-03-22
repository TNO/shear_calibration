"""
    Testing the python package: pacal.
"""

import pacal as p

p.params.general.parallel = False


G = p.GumbelDistr(mu=0, sigma=1)
N1 = p.NormalDistr(mu=0.0, sigma=1.0)
LN1 = p.exp(N1)
N2 = p.NormalDistr(mu=1.0, sigma=0.1)

T = G * LN1 * N2

T.plot()
T.summary()
