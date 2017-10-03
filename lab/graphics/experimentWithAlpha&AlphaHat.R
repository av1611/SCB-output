rm (list = ls())
myPath <- "out"
if(!dir.exists(myPath))
{
  dir.create(myPath)
}
myFileName <- gsub(" ","_",paste("myGraphWithAlpha&AlphaHat_",Sys.time(),".jpg",sep = ""))
myFileName <- gsub(":","_",myFileName)
#jpeg(paste(myPath,"/",myFileName,sep=""))

alpha1 <- 0.8
alphaHat1 <- c(0.72,0.75,0.82,0.78,0.7,0.83)
alpha2 <- 0.6
alphaHat2 <- c(0.57,0.58,0.63,0.61,0.6,0.59)
alpha3 <- 0.4
alphaHat3 <- c(0.38,0.39,0.43,0.43,0.4,0.39)
alpha4 <- 0.2
alphaHat4 <- c(0.18,0.19,0.23,0.27,0.17,0.2)

arrayOfXAplphaHat1 <- array(alpha1,dim = (length(alphaHat1)))
arrayOfXAplphaHat2 <- array(alpha2,dim = (length(alphaHat2)))
arrayOfXAplphaHat3 <- array(alpha3,dim = (length(alphaHat3)))
arrayOfXAplphaHat4 <- array(alpha4,dim = (length(alphaHat4)))

plot(x=c(0,1),y=c(0,1),type ="n",main="Alpha and AlphaHat")
abline(a = 0 , b =1)

points(x=alpha1, y=0,col="blue")
points(x=arrayOfXAplphaHat1,y=alphaHat1,col="red")

points(x=alpha2, y=0,col="blue")
points(x=arrayOfXAplphaHat2,y=alphaHat2,col="red")

points(x=alpha3, y=0,col="blue")
points(x=arrayOfXAplphaHat3,y=alphaHat3,col="red")

points(x=alpha4, y=0,col="blue")
points(x=arrayOfXAplphaHat4,y=alphaHat4,col="red")

legend("topleft",c("Alpha","AplhaHat"),fill = c("blue","red"))

#dev.off()
