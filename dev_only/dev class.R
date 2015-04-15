library(MALDIquant)

# object oriented approach for tof measurements

setClass("TofMeasurement",
         slots=list(metaDataScan="list", file="character",
                    .indexHelp="list", .fid="H5IdComponent", .tofBlock="H5IdComponent"),
         prototype=list(metaDataScan=list(), file="",
                        .indexHelp=list())
         )

# constructor
createTofMeasurement <- function(file){
  .fid <-H5Fopen(tof.h5)
  scan.attr <- h5readAttributes(.fid,'FullSpectra')
  
  .tofBlock <- get.raw.tofblock(.fid)
  h5spaceFile <- H5Dget_space(tofblock)
  dims <- H5Sget_simple_extent_dims(h5spaceFile)
  calc.indices <- expand.grid(buf=1:dims$size[[3]], write=1:dims$size[[4]])
  h5spaceMem <- H5Screate_simple(dims$size[[1]])
  .indexHelp <- list(N=nrow(calc.indices),
       h5spaceFile=h5spaceFile, 
       h5spaceMem=h5spaceMem, 
       calc.indices=calc.indices, 
       dims=dims$size)
  
  new(Class = "TofMeasurement", file=file, metaDataScan=scan.attr, 
      .fid=.fid, .indexHelp = .indexHelp, .tofBlock=.tofBlock)
}

setMethod(f="getScan", 
          signature = signature("TofMeasurement"),
          definition=function(object, index=1){
            read.spec.ind(object@.tofBlock, object@.indexHelp,  index) 
          })
