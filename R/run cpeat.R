
cpeat<- read.csv("cpeat.csv")


nls <- clymo(cpeat, 0.01, 0.01)

cpeat <- prep(cpeat, nls)

acro <- acrotelm(cpeat)

mega <- mega_bog(cpeat)



