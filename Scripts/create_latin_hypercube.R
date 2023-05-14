#### Creating a latin hypercube sampling design

# we run a combination of densities, swimming speed and chemodetection ranges
library(DoE.wrapper)


Hunting.LHS <- lhs.design(seed = 42,
  nruns = 100,
  nfactors = 2,
  type = "maximin",
  digits = c(8,8),
  factor.names = list(
    "glaucus_Chemodetection" = c(0, 0.05),
    "glaucus_Speed" = c(0, 0.05)))

write_csv(Hunting.LHS, './Outputs/LHS_Hunting_Behaviour.csv')
