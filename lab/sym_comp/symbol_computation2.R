rm (list = ls())

if (! require ("Ryacas"))
  install.packages("Ryacas")
library(Ryacas)

x <- Sym("x")

phi <- function(x)
{
  (1 / sqrt(2 * pi)) * exp((-x ^ 2) / 2)
}

normDifKernel <- function(x)
{
  2 * phi(x) - phi(x / sqrt(2)) * (1 / sqrt(2))
}

kernelDeriv <- normDifKernel
body(kernelDeriv) <- yacas((deriv(normDifKernel(x))))[[1]]

kernelDerivSquared <- normDifKernel
body(kernelDerivSquared) <- yacas(kernelDeriv(x) ^ 2)[[1]]

lowBound <- -6
upBound <- 6

intSqDer <- integrate(f = kernelDerivSquared,
                      lower = lowBound,
                      upper = upBound)[[1]]
intSqDer


kernelSquared <- normDifKernel
body(kernelSquared) <- yacas(normDifKernel(x) ^ 2)[[1]]

phiKNormalDif <- integrate(f = kernelSquared,
                           lower = lowBound,
                           upper = upBound)[[1]]
phiKNormalDif


CK <- 0.5 * log (intSqDer / (4 * pi ^ 2 * phiKNormalDif))
CK
