### load the package and the data for the BCG vaccine meta-analysis
#install.packages('readr')
library("metafor")

### calculate the log relative risks and corresponding sampling variances
data <- read.csv('data.csv')
print(data)
dat <- escalc(measure = "OR", ai = tpos, bi = tneg, ci = cpos,
              di = cneg, data = data, append = TRUE)
print(dat[,-c(4:7)], row.names = FALSE)


k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
dat.fm

print('--- Random effect model ---')
res <- rma(ai = tpos, bi = tneg, ci = cpos,
           di = cneg, data = dat, measure = "OR")
res

### confidence intervals for tau^2, I^2, and H^2

confint(res)

print('---------------sensitivity')
one_out = leave1out(res)
print(one_out)


### forst plot

forest(res, slab = paste(dat$author, dat$year, sep = ", "),
       xlim = c(-16, 6), at = log(c(.05, .25, 1, 4)), atransf = exp,
       ilab = cbind(dat$tpos, dat$tneg, dat$cpos, dat$cneg),
       ilab.xpos = c(-9.5, -8, -6, -4.5), cex = .75)
op <- par(cex = .75, font = 2)
text(c(-9.5, -8, -6, -4.5), 15, c("HPV+", "HPV-", "HPV+", "HPV-"))
text(c(-8.75, -5.25),       16, c("OLP", "Control"))
text(-16,                   15, "Author(s) and Year",     pos = 4)
text(6,                     15, "OR [95% CI]", pos = 2)
par(op)

###########3 Egger's test
('-----------egger test fixed')
regtest(x=res, vi, sei, data=dat, model="lm", predictor="sei", ret.fit=FALSE)
('-----------egger test random')
regtest(x=res, vi, sei, data=dat, model="rma", predictor="sei", ret.fit=FALSE)
('-----------egger test- return fitted')
regtest(x=res, vi, sei, data=dat, model="rma", predictor="sei", ret.fit=TRUE)

print('_________ funnel plots')
res <- rma(yi, vi, data = dat)
funnel(res, main = "Random-Effects Model")

### trim and fill method
print('----random with fill and trim----')

res <- rma(ai = tpos, bi = tneg, ci = cpos,
           di = cneg, data = dat, measure = "OR")
rtf <- trimfill(res)
print('--- Trimmed and filled result ---')
rtf
print('--- Funnel plot for filled and trimmed ---')
funnel(rtf, main = "Random-Effects Model, Fill-and-trim")

filled <- data.frame(yi = rtf$yi, vi = rtf$vi, fill = rtf$fill)
filled


forest(rtf, xlim = c(-16, 14), at = log(c(.005, .05, .25, 1, 10, 100, 1000)), atransf = exp, cex = .75,
       slab = append(paste(dat$author, dat$year, sep = ", "),
       c("Filled 1", "Filled 2", "Filled 3", "Filled 4", "Filled 5", "Filled 6")))
