install.packages("Derive")
library(Derive)

x <- Sym("x")

K <- function(x) {
  x^3
}

KDerive <- Ryacas::yacas(deriv(K(x), x))
KDerive$text

KDeriveAsChExpr <- Ryacas::as.character.Expr(KDerive$text)
# KDeriveAsChExpr <- paste0("Integrate(x)", KDeriveAsChExpr)


lowBound <- -1
upBound <- 1

KDerivIntegrated <- Ryacas::Integrate.default(f = KDerive$text,
                                              x = x,
                                              a = lowBound,
                                              b = upBound)



print(yacas(KDerivIntegrated))


