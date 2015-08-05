library(rhdf5)

setClass("TofH5", 
         representation(file.name='character',
         				n.scans='numeric', 
         				meas.reader='H5IdComponent',
         				.tofbloc='H5IdComponent',
         				.indexhelp='list',
         				.Data='list') 
         )
 
