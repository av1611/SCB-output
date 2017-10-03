

if (! require ("Ryacas"))
  install.packages("Ryacas")
library(Ryacas)


K = function (x)
{
  return (2*(2*pi)^(-1/2) * exp (-(x^2) / 2) - (2*pi)^(-1/2) *
                            exp (-((x/sqrt(2))^2) / 2) / sqrt (2))
}



INT_SQ_DER = 0.306951
PHI_K_NORMAL_DIF = 0.4065

CK = 0.5 * log (INT_SQ_DER / (4 * pi^2 * PHI_K_NORMAL_DIF))
CK


# AVG: Taking a derivative in a sym computation via Ryacas
if (! require ("Ryacas"))
  install.packages("Ryacas")
library(Ryacas)

x <- Sym("x")

K = function (x) {
  2*(2*pi)^(-1/2) * exp (-(x^2) / 2) -
    (2*pi)^(-1/2) * exp (-((x/sqrt(2))^2) / 2) / sqrt (2)
}


KDerive <- Ryacas::yacas(deriv(K(x), x))
KDerive$text
# expression(-0.7978845608 * (exp(-x^2/2) * (2 * x))/2 - -0.3989422804 *
#              (exp(-(x/1.4142135623731)^2/2) * (2 * x))/5.6568542495)

lowBound <- -1
upBound <- 1

KDerivIntegrated <- Ryacas::Integrate.default(f = KDerive$text,
                                              x = x,
                                              a = lowBound,
                                              b = upBound)

as.character.Sym(KDerivIntegrated)



kDerivIndefiniteIntegral <- Ryacas::Integrate.default(f = KDerive$text,
                                                      x = x)

print(yacas(kDerivIndefiniteIntegral))
