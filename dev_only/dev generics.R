# create generics

if (is.null(getGeneric("length"))){
  setGeneric("length", function(x) standardGeneric("length"))
}

# 
if (is.null(getGeneric("estimateNoise"))) {
  setGeneric("estimateNoise",
             function(object, ...) standardGeneric("estimateNoise"))
}


if (is.null(getGeneric("getScan"))) {
  setGeneric("getScan",
             function(object, index) standardGeneric("getScan"))
}
