# object oriented approach for tof measurements

#' Class TofMeasurement
#' 
#' Class \code{TofMeasurement} handles access to a Tof hdf5 file
#' 
#' @slot file Name of the hdf5 file
#' @slot metaDataScan some extra informations
#' @import rhdf5
#' @name tofmeasurement
#' @exportClass TofMeasurement
TofMeasurement <- setClass("TofMeasurement",
         slots=list(metaDataScan="list", file="character",
                    .indexHelp="list", .fid="H5IdComponent", .tofBlock="H5IdComponent"),
         prototype=list(metaDataScan=list(), file="",
                        .indexHelp=list())
         )

#' constructor for TofMeasurement objects
#' 
#'  @param file the name of the hdf5 file containing the scan
#'  
#'  @return a TofMeasurement object representing the file
#'  
#'  @usage createTofMeasurement(file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5"))
#'  
#'  @export
#'  @docType methods
createTofMeasurement <- function(file){
  .fid <- H5Fopen(tof.h5)
  scan.attr <- h5readAttributes(.fid,'FullSpectra')
  
  .tofBlock <- get.raw.tofblock(.fid)
  h5spaceFile <- H5Dget_space(.tofBlock)
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

