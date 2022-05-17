library(tidyverse)

## sigmoid decay curve
# 
# y = (1 / (1 + exp((b * x / c) - a))) * w
# 
# a - shifts curve to right or left 
# b - determines spread of curve, or slope of the rapidly decreasing part of
#     curve
# c - scalar to adjust total distance of interest (=distance in meters divided
#     by 20)
# x - distance in meters from threat 
# w - weight of threat (maximum value) 
# a:b ratio determines the inflection point
# 10:1 ratio for inflection at 500 m

## example
a <- 5
b <- 0.5
c <- 50
w <- 1000
x <- c(1:1000)
val <- (1 / (1 + exp((b * x / c) - a))) * w

plot(x,val)
exp(1)

## generate a sigmoid curve from a list of parameters
makeCurve <- function(params) {
  points <- (1 / 
               (1 + exp((params[2] * c(1:1000) / 
                           params[3]) - params[1]))) * params[4]
  points
}


## Steepness: increase a and b in tandem to steepen the curve while maintaining
## inflection
## make curves
params1 <- c(1, 0.1, 50, 1000)
curve1 <- makeCurve(params1)
params2 <- c(3, 0.3, 50, 1000)
curve2 <- makeCurve(params2)
params3 <- c(5, 0.5, 50, 1000)
curve3 <- makeCurve(params3)
params4 <- c(7, 0.7, 50, 1000)
curve4 <- makeCurve(params4)
params5 <- c(9, 0.9, 50, 1000)
curve5 <- makeCurve(params5)

## plot curves
ggplot() + 
  geom_point(aes(x = c(1:1000), y = curve1), color = 'red') +
  geom_point(aes(x = c(1:1000), y = curve2), color = 'orange') +
  geom_point(aes(x = c(1:1000), y = curve3), color = 'green') +
  geom_point(aes(x = c(1:1000), y = curve4), color = 'blue') +
  geom_point(aes(x = c(1:1000), y = curve5), color = 'purple') +
  labs(y = 'Weight', x = 'Distance (m)')
  
## Inflection: increase a to shift inflection farther away
## make curves
params1 <- c(1, 0.5, 50, 1000)
curve1 <- makeCurve(params1)
params2 <- c(3, 0.5, 50, 1000)
curve2 <- makeCurve(params2)
params3 <- c(5, 0.5, 50, 1000)
curve3 <- makeCurve(params3)
params4 <- c(7, 0.5, 50, 1000)
curve4 <- makeCurve(params4)
params5 <- c(9, 0.5, 50, 1000)
curve5 <- makeCurve(params5)

## plot curves
ggplot() + 
  geom_point(aes(x = c(1:1000), y = curve1), color = 'red') +
  geom_point(aes(x = c(1:1000), y = curve2), color = 'orange') +
  geom_point(aes(x = c(1:1000), y = curve3), color = 'green') +
  geom_point(aes(x = c(1:1000), y = curve4), color = 'blue') +
  geom_point(aes(x = c(1:1000), y = curve5), color = 'purple') +
  labs(y = 'Weight', x = 'Distance (m)')

## Stretch: increase b to decrease length of decay
## make curves
params1 <- c(5, 0.1, 50, 1000)
curve1 <- makeCurve(params1)
params2 <- c(5, 0.3, 50, 1000)
curve2 <- makeCurve(params2)
params3 <- c(5, 0.5, 50, 1000)
curve3 <- makeCurve(params3)
params4 <- c(5, 0.7, 50, 1000)
curve4 <- makeCurve(params4)
params5 <- c(5, 0.9, 50, 1000)
curve5 <- makeCurve(params5)

## plot curves
ggplot() + 
  geom_point(aes(x = c(1:1000), y = curve1), color = 'red') +
  geom_point(aes(x = c(1:1000), y = curve2), color = 'orange') +
  geom_point(aes(x = c(1:1000), y = curve3), color = 'green') +
  geom_point(aes(x = c(1:1000), y = curve4), color = 'blue') +
  geom_point(aes(x = c(1:1000), y = curve5), color = 'purple') +
  labs(y = 'Weight', x = 'Distance (m)')

## params for initial trials - future trials could look into lit for
## justification? How else could we fit these? Examine ring coeffs?
## make curves
params1 <- c(4.5, 3.5, 50, 1000) # 1.2857
curve1 <- makeCurve(params1)
params2 <- c(6, 2, 50, 1000) # 3
curve2 <- makeCurve(params2)
params3 <- c(8, 1.25, 50, 1000) # 6.4
curve3 <- makeCurve(params3)
params4 <- c(10, 1, 50, 1000) # 1
curve4 <- makeCurve(params4)
params5 <- c(11, 0.75, 50, 1000) # 14.6666
curve5 <- makeCurve(params5)

## plot curves
ggplot() + 
  geom_point(aes(x = c(1:1000), y = curve1), color = 'red') +
  geom_point(aes(x = c(1:1000), y = curve2), color = 'orange') +
  geom_point(aes(x = c(1:1000), y = curve3), color = 'green') +
  geom_point(aes(x = c(1:1000), y = curve4), color = 'blue') +
  geom_point(aes(x = c(1:1000), y = curve5), color = 'purple') +
  labs(y = 'Weight', x = 'Distance (m)')
