
############ First.lib ###############

.onLoad <- function(lib, pkg){
   require(haplo.stats)
   require(survival)
   require(mvtnorm)
   library.dynam("SNPassoc", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("SNPassoc", libpath)


############ End of .First.lib ###############


