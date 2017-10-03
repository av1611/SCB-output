# Define K() = x^3
# Differentiate it =  K'
# Square it =  (K')^2
# simplify
# Integrate it from -5 to 5
# Tex it
# Export to file
# (x^2)^3 = x^6

if (! require ("Ryacas"))
  install.packages("Ryacas")
library(Ryacas)

x <- Sym("x")

K <- function(x) {
  x^3
}

KDerive <- Ryacas::yacas(deriv(K(x), x))
KDerive$text
lowBound <- -1
upBound <- 1

KDerivIntegrated <- Ryacas::Integrate.default(f = KDerive$text,
                                              x = x,
                                              a = lowBound,
                                              b = upBound)

as.numeric (KDerivIntegrated)
as.character.Sym (KDerivIntegrated)

print(yacas(KDerivIntegrated))

kDerivIndefiniteIntegral <- Ryacas::Integrate.default(f = KDerive$tex,
                                              x = x
                                              )

print(yacas(kDerivIndefiniteIntegral))

kDerivIndefiniteIntegralAsCh <- Ryacas::as.character.Expr(kDerivIndefiniteIntegral)
kDerivIndefiniteIntegralAsCh

yacas("myExprCubeDerivAsCh := 3 * x^2")
yacas("TeXForm(myExprCubeDerivAsCh)", retclass = "unquote")


