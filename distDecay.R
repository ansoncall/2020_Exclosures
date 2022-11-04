# sigmoid decay funcs
dist <- seq(1,1000)
sig1 <- (1.0 / (1.0 + 2.718282 ^ ((3.5 * dist / 50) - 4.5)))
sig1 <- (1.0 / (1.0 + 2.718282 ^ ((3.5 * dist / 50) - 4.5)))
sig2 <- 1/((1+2.718)^((2.0*dist/50)-6.0))
sig3 <- 1/((1+2.718)^((1.25*dist/50)-8.0))

df <- data.frame(dist, sig1, sig2, sig3)

ggplot(df) +
  geom_point(aes(dist, sig1), color = 'red')#+
  geom_point(aes(dist, sig2), color = 'blue')+
  geom_point(aes(dist, sig3), color = 'orange')
