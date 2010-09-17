".First.lib" <-function(lib, pkg)
{
  library.dynam("baybi", package = pkg, lib.loc = lib)
  return(invisible(0)) 
}
