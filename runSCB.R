rm (list = ls())

if (!require("rstudioapi"))
  install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (!require("roxygen2"))
  install.packages("roxygen2")
library(roxygen2)

if (!require("devtools"))
  install.packages("devtools")
library(devtools)

if (!require("ggplot2"))
  install.packages("ggplot2")
require(ggplot2)

if (!require("gridExtra"))
  install.packages("gridExtra")
library(devtools)
require(gridExtra)

# if (!require("KernSmooth"))
#   install.packages("KernSmooth")
# library(KernSmooth)
# require(KernSmooth)

setwd("SCB/") # error? at this level there is no more scb

roxygen2::roxygenise()
devtools::document()

# Loading the functions


noiseMean = 0
noiseSd = 1
sampleSize = 100
corruption = sin #function(x){sqrt(x)}
T <- createT(sampleSize = sampleSize)
T
Psi <- createPsi(model = corruption,
                 time = T)
Psi
Z = createZ(sampleSize = sampleSize,
            mean = noiseMean,
            sd = noiseSd)
Z

X <- createX(model = ma1,
             sampleSize = sampleSize,
             Z = Z,
             psi = Psi)
X
corruptedMA1 = data.frame(T*sampleSize,X)
colnames(corruptedMA1) <- c('time','sample')

pl_title = sprintf("Noise distribution N(%d,%.2f), corruption is sin, time steps %d", noiseMean, noiseSd, sampleSize)


noise_mean = 0
noise_sd = 1
sampling_time = 350
corruption = sin #function(x){sqrt(x)}
noise <- createSleeper(number_steps = sampling_time)
noise

psi <- createPsi(model = corruption,
                 sleeper = sleeper)
psi

noise = createNoise(time = sampling_time,
            mean = noise_mean,
            sd = noise_sd)
noise

sample <- createSample(model = ma1,
             time = sampling_time,
             noise = noise,
             psi = psi)
sample

corrupted_ma1 = data.frame(T*sampling_time, sample)
colnames(corrupted_ma1) <- c('time','sample')

pl_title = sprintf("Noise distribution N(%d,%.2f), corruption is sin, time steps %d", noise_mean, noise_sd, sampling_time)
if (! dir.exists ("Plots"))
  dir.create("Plots")
setwd ("Plots")

jpeg(filename = paste0(pl_title, ".jpg"), width = 2*585, height = 2*403, units = "px", quality = 80)
ggplot(data=corruptedMA1, aes(x = time, y = sample)) + geom_line() +
  xlab("Time steps") +
  ylab("Sample path") +
  ggtitle(pl_title) +
  theme(plot.title = element_text(face = "bold", hjust = .5, size = 25), axis.text = element_text(size=20), axis.title = element_text(size = 20))
graphics.off()
# dev.off()
setwd ("..")
getwd()


jpeg (filename = "tmp0.jpeg")
dev.off()
<<<<<<< HEAD

# W-process generation & plotting
gaussKernel <- function(x){
  1/sqrt(2*pi)*exp(-0.5*x^2)
}
sigmoidKernel <- function(x){
  2/pi * 1/(exp(x) + exp(-x))
}
W <- data.frame(cbind(1:sampleSize, createW(kernel=sigmoidKernel,bandwidth=0.5,sampleSize)))
colnames(W) <- c('time', 'W')
ggplot(data = W, aes(x=time, y=W)) + geom_line() + ggtitle('W kernel process') +
  theme(plot.title = element_text(face = "bold", hjust = .5, size = 20), axis.text = element_text(size=15), axis.title = element_text(size = 15))

# Bootstraps for MA1
psiFixed = 0.5
bSampleSize = 10
bandwidth = 0.5
sourceMA1 = ma1(sampleSize, psi = rep(psiFixed, sampleSize), rnorm(sampleSize,0,1))
bootstrappedSamples = bootstrapOnMA1(sourceMA1, bSampleSize = bSampleSize, kernel = sigmoidKernel, bandwidth = bandwidth)

b1 = data.frame(cbind(1:sampleSize, bootstrappedSamples[1, ]))
colnames(b1) <- c('time', 'sample')

b2 = data.frame(cbind(1:sampleSize, bootstrappedSamples[2, ]))
colnames(b2) <- c('time', 'sample')

bootstr1 = ggplot(b1, aes(x=time, y=sample)) + geom_line()
bootstr2 = ggplot(b2, aes(x=time, y=sample)) + geom_line()
grid.arrange(bootstr1, bootstr2, ncol=2, nrow=1)

=======
gauss_ker <- function(x){
  1/sqrt(2*pi)*exp(-0.5*x^2)
}
W <- data.frame(cbind(1:sampling_time, createBootstrapMultiplier(kernel=gauss_ker, bandwidth=0.5, sampling_time)))
colnames(W) <- c('time', 'W')
ggplot(data = W, aes(x=time, y=W)) + geom_line() + ggtitle('W kernel process') +
  theme(plot.title = element_text(face = "bold", hjust = .5, size = 20), axis.text = element_text(size=15), axis.title = element_text(size = 15))
>>>>>>> 4bfb282544601b59851f24cef9e00599c5b92b58

