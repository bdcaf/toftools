#" open TOF file
#"
#" open required for reading tof data
#" @export
#" @docType methods
#" @rdname open-methods
#" @exportMethod open
setGeneric("open", function(object) {
  standardGeneric("open")
})

#" @rdname open-methods
setMethod("open", signature(object = "tof_h5"), function(object) {
                        object@.meas.reader <- H5Fopen(object@file.name)
                        object@.tofblock <- raw_tofblock(object@.meas.reader)
                        object@.indexhelp <- tof.indexhelp(object@.tofblock)
      object
})


#" close TOF file
#"
#" close usually not required
#" @export
#" @docType methods
#" @rdname close-methods
#" @exportMethod close
setGeneric("close", function(object) {
  standardGeneric("close")
})

#" @rdname close-methods
setMethod("close", signature(object = "tof_h5"), function(object) {
                        H5Fclose(object@.meas.reader)

                        object
})


