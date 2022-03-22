"""
    Testing the python package: pacal.
"""

import pacal as p

p.params.general.parallel = False

N = p.NormalDistr()
N.plot()
p.show()

N.cdf(0)  # should produce 0.5
N.summary()  # prints summary: mean, variance, etc.

C = N / N  # Cauchy distribution
C.summary()

L = N * N + N * N  # Laplace distribution
L.plot()
p.show()


G = p.GumbelDistr(mu=0, sigma=1)
N = p.NormalDistr(mu=0.0, sigma=1.0)
LN = p.exp(N)

T = G * LN

T.plot()

T.pdf(1)


# Define a new distribution type

# FunDistr(f, breakPoints = [0, Inf])
