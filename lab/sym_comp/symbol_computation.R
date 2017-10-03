if (! require ("Ryacas"))
install.packages("Ryacas")
library(Ryacas)

x <- Sym("x")

myExprCubeX <- function(x) {
  x^3
}

myExpr <- function(x) {
  returnValue <- 2 *
    ((2 * pi) ^ (-1/2) *
       exp(-(x^2) / 2) -  (2*pi)^(-1/2) *
       exp (-((x/sqrt(2))^2) / 2) / sqrt(2))^2
}






myExprCubeDeriv <- Ryacas::yacas(deriv(myExprCubeX(x), x))
myExprCubeDeriv$text

yacas("myExprCubeDerivAsCh := 3 * x^2")
yacas("TeXForm(myExprCubeDerivAsCh)", retclass = "unquote")

# $3 x ^{2}$

yacas("PrettyForm(myExprCubeDerivAsCh)")



myExprDeriv <- Ryacas::yacas(deriv(myExpr(x), x))
myExprDeriv$text

yacas("myExprDerivAsCh := 4 * ((-0.3989422804 * (exp(-x^2/2) * (2 * x))/2 -
    -0.3989422804 * (exp(-(x/1.4142135623731)^2/2) * (2 * x))/5.6568542495) *
    (0.398942280401433 * exp(-x^2/2) - 0.398942280401433 * exp(-(x/1.4142135623731)^2/2)/1.4142135623731))")
yacas("TeXForm(myExprDerivAsCh)", retclass = "unquote")

# $4 \left( \frac{0.3989422804 \mathrm{ exp }\left(  - \frac{\left( \frac{x}{1.4142135623731} \right)  ^{2}}{2} \right)  2 x}{5.6568542495}  - \frac{0.3989422804 \mathrm{ exp }\left(  - \frac{x ^{2}}{2} \right)  2 x}{2} \right)  \left( 0.398942280401433 \mathrm{ exp }\left(  - \frac{x ^{2}}{2} \right)  - \frac{0.398942280401433 \mathrm{ exp }\left(  - \frac{\left( \frac{x}{1.4142135623731} \right)  ^{2}}{2} \right) }{1.4142135623731} \right) $


yacas("PrettyForm(myExprDerivAsCh)")



