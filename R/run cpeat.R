cpeat_data<- read.csv("cpeat.csv")

library(usethis)


nls <- clymo(cpeat, 0.01, 0.01)

cpeat <- prep(cpeat, nls)

acro <- acrotelm(cpeat)

mega <- mega_bog(cpeat)


use_data(cpeat_data)
