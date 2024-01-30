getwd()
### TIME SERIES ANALYSES
library(forecast)
library(tseries)

### CROSS-WAVELET ANALYSES
library(WaveletComp)
library(dplyr)
library(matrixStats)
library(tidyr)
library(ggplot2)

### Example with synthetic data
# The following example is modified from Veleda et al, 2012:
series.length <- 3*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 2*24, length = series.length)
x3 <- periodic.series(start.period = 4*24, length = series.length)
x4 <- periodic.series(start.period = 8*24, length = series.length)
x5 <- periodic.series(start.period = 16*24, length = series.length)
x6 <- periodic.series(start.period = 32*24, length = series.length)
x7 <- periodic.series(start.period = 64*24, length = series.length)
x8 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3 + x4 + 3*x5 + x6 + x7 + x8 + rnorm(series.length)
y <- x1 + x2 + x3 + x4 - 3*x5 + x6 + 3*x7 + x8 + rnorm(series.length)

matplot(data.frame(x, y), type = "l", lty = 1, xaxs = "i", col = 1:2, 
        xlab = "index", ylab = "",
        main = "hourly series with periods of 1, 2, 4, 8, 16, 32, 64, 128 days", 
        sub = "(out of phase at period 16, different amplitudes at period 64)")
legend("topright", legend = c("x","y"), col = 1:2, lty = 1)

## The following dates refer to the local time zone 
## (possibly allowing for daylight saving time):      
my.date <- seq(as.POSIXct("2014-10-14 00:00:00", format = "%F %T"), 
               by = "hour", 
               length.out = series.length)     
my.data <- data.frame(date = my.date, x = x, y = y)

## Computation of cross-wavelet power and wavelet coherence, x over y:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24',
## resulting in one time unit equaling one day. 
## This is also the time unit in which periods are measured.
## There is an option to store the date format and time zone as additional
## parameters within object 'my.wc' for later reference. 

my.wc <- analyze.coherency(my.data, c("x","y"), 
                           loess.span = 0, 
                           dt = 1/24, dj = 1/20, 
                           window.size.t = 1, window.size.s = 1/2, 
                           lowerPeriod = 1/4,
                           make.pval = TRUE, n.sim = 10,
                           date.format = "%F %T", date.tz = "")
## Note:                           
## By default, Bartlett windows are used for smoothing in order to obtain
## the wavelet coherence spectrum; window lengths in this example:
## 1*24 + 1 = 25 observations in time direction,
## (1/2)*20 + 1 = 11 units in scale (period) direction.                             

## Plot of cross-wavelet power 
## (with color breakpoints according to quantiles):
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)")

## The same plot, now with calendar axis
## (according to date format stored in 'my.wc'):
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)", show.date = TRUE)   

## Plot of average cross-wavelet power:
wc.avg(my.wc, siglvl = 0.05, sigcol = 'red', 
       periodlab = "period (days)")

## Plot of wavelet coherence 
## (with color breakpoints according to quantiles):
wc.image(my.wc, which.image = "wc",  main = "wavelet coherence, x over y", 
         legend.params = list(lab = "wavelet coherence levels", 
                              lab.line = 3.5, label.digits = 3),
         periodlab = "period (days)")

## plot of average coherence:
wc.avg(my.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (days)")



### Now with my data
## Load data
# Load surface water level time series (daily)
df_srs2.level <- read.csv("fce.srs2.level.csv")
plot.ts(df_srs2.level$y)

# Load surface water dissolved organic carbon time series (monthly)
df_srs2.doc <-read.csv("fce.srs2.doc.csv")
plot.ts(df_srs2.doc$y)

# See WaveletComp tour at: http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf
## wavelets and cross-wavelet
# wavelet on flow data
?analyze.wavelet
head(df_srs2.level)
my.wt <- analyze.wavelet(df_srs2.level, "y",make.pval = TRUE, n.sim = 2) # n.simulations will need to be higher for the final figure, but takes very long
wt.image(my.wt, main = "surface water level",
         periodlab = "period (months)",
         label.time.axis = T, show.date = T, date.format = "%m/%d/%y",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))


# wavelet on surface water dissolved organic carbon
my.wt <- analyze.wavelet(df_srs2.doc, "y",make.pval = TRUE, n.sim = 2) # n.simulations will need to be higher for the final figure, but takes very long
wt.image(my.wt, main = "surface water DOC",
         periodlab = "period (months)",
         label.time.axis = T, show.date = T, date.format = "%m/%d/%y",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))


# cross-wavelet flow & temp (contours show "shared" scales of environmental variation)
my.data<-cbind(df_srs2.level, df_srs2.doc$y)# merge data
colnames(my.data)<-c("date","level","doc")
head(my.data)
my.wc <- analyze.coherency(my.data, c("level","doc"), n.sim = 2) # n.simulations will need to be higher for the final figure, but takes very long
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (months)")

# plot of average coherence:
wc.avg(my.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (monthss)")
