
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("SNPassoc", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("SNPassoc", libpath)


############ End of .First.lib ###############


