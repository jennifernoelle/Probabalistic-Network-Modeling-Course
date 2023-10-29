library(mcclust.ext)
data(galaxy.fit)
x=data.frame(x=galaxy.fit$x)
data(galaxy.pred)
data(galaxy.draw)
# Find representative partition of posterior
# Variation of Information (minimizes lower bound to VI)
psm=comp.psm(galaxy.draw)
galaxy.VI=minVI(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.VI)
plot(galaxy.VI,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Compute Variation of Information
VI(galaxy.VI$cl,galaxy.draw)
# Binder
galaxy.B=minbinder.ext(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.B)
plot(galaxy.B,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Uncertainty in partition estimate
galaxy.cb=credibleball(galaxy.VI$cl[1,],galaxy.draw)
summary(galaxy.cb)
plot(galaxy.cb,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Compare with uncertainty in heat map of posterior similarity matrix
plotpsm(psm)

p = seq(0, 10, length=100)

# Create plot of Beta distribution with various shape parameters
plot(p, dgamma(p, 1, 1), type='l', ylab = 'Density')
lines(p, dgamma(p, 1.25, 1), col='red') 
lines(p, dgamma(p, 1.5, 1), col='blue')

legend(x = .7, y = 4, legend = c('Beta(1, 10)','Beta(2, 2)','Beta(1,1)'),
       lty=c(1,1,1), col=c('black', 'red', 'blue'))
