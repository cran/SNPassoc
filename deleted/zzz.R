
############ First.lib ###############

.First.lib <- function(lib, pkg){
   require(haplo.stats)
   require(survival)
   library.dynam("SNPassoc", pkg, lib)
}
############ End of .First.lib ###############


