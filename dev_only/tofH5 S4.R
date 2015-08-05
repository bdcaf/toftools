library(rhdf5)

setClass("TofH5", 
         representation(file.name='character',
         				n.scans='numeric', 
         				meas.reader='H5IdComponent',
         				.Data='list') 
         )
 
# open tof file
setGeneric("open", function(object) {
  standardGeneric("open")
})


setMethod('open', signature(object = 'TofH5'), function(object) {
			object@meas.reader <- H5Fopen(object@file.name)
})


