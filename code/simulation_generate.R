# devtools::install_github('xzhoulab/SRTsim')
library(SRTsim)

setwd("../data/sim_data")
example_count   <- exampleLIBD$count
example_loc     <- exampleLIBD$info[,c("imagecol","imagerow","layer")]
colnames(example_loc) <- c("x","y","label")
simSRT  <- createSRT(count_in=example_count,loc_in =example_loc)


set.seed(1)
simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
simSRT1 <- srtsim_count(simSRT1)

set.seed(1)
simSRT2 <- srtsim_fit(simSRT,sim_schem="domain")
simSRT2 <- srtsim_count(simSRT2)

write.csv(as.matrix(simSRT1@simCounts), "ref_based/ref_based_tissue_X.csv")
write.csv(as.matrix(simSRT2@simCounts), "ref_based/ref_based_domain_X.csv")
location1 = data.frame("x" = simSRT1@simcolData@listData[[1]], "y" = simSRT1@simcolData@listData[[2]], "label" = simSRT1@simcolData@listData[[3]])
location2 = data.frame("x" = simSRT2@simcolData@listData[[1]], "y" = simSRT2@simcolData@listData[[2]], "label" = simSRT2@simcolData@listData[[3]])
write.csv(location1, "ref_based/ref_based_tissue_S.csv")
write.csv(location2, "ref_based/ref_based_domain_S.csv")






#### reference-free simulation 
shinySRT1 = SRTsim_shiny()
# number of points: 
# number of clusters: 2965
# mean fold: 1,2,3,5
# number of higher gene: 100
# number of lower gene: 500
# number of noise gene: 1400
# zero proportion: 0.1
# dispersion: 0.1
# mean: 1
setwd("../data/sim_data")
save(shinySRT1, file = "ref_free/shiny.rdata")
load("ref_free/shiny.rdata")
write.csv(shinySRT1$simInfo, file = "ref_free/S.csv")

rep = 20
dispersion = c(0.1, 0.5, 0.9, 1.3, 1.8)

for(i in 1:length(dispersion))
{
  tmp = shinySRT1
  tmp$simcountParam$Dispersion = dispersion[i]
  for(j in 1:rep)
  {
    count_new = reGenCountshiny(tmp, NewSeed = 10*j)
    filepath = paste0("ref_free/X_disp_", dispersion[i], "_rep_", j, ".csv")
    write.csv(count_new, file = filepath)
  }
}

rep = shinySRT1
rep$simcountParam$mu = 4
rep1 = reGenCountshiny(rep)










