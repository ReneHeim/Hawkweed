#' Remove spectral outliers via depth measure
#'
#' @param specdf A df containing spectral data with the first col named "Type"
#' and containing the categorical response variables and all following columns
#' named according to the spectral band number (e.g. 450, 455, 1640...)
#' @dfuncstr A string indicating the depth measure function
#' @nbno Number of bootstrap samples
#' @smono Number indicating smooting factor
#' @trimno Number indicating trim value
#' @return The same df as provided to the function with all outliers removed
#' @examples
#' rmv.func.outlier(specdf)
#' 


rmv.funct.outlier <- function(specdf, dfuncstr, nbno, smono, trimno){
        
        require(fda)
        require(fda.usc)
  
  dat.mat <-
    as.matrix(specdf[, 2:length(names(specdf))]) # As matrix to be able to transform the object
  
  labnames <- list(main="Spectra", 
                   xlab="Wavelength [nm]", 
                   ylab="Reflectance [%]")
  
  myfdata <-
    fdata(dat.mat,
          argvals = as.integer(names(specdf[, 2:length(names(specdf))])),
          names = labnames) # fda pkg needs fdata object to run, created here.
  
  outlier.mat <-
    outliers.depth.trim(
      myfdata,
      dfunc = dfuncstr,
      nb = nbno,
      smo = smono,
      trim = trimno
    )
  
  outlier.vector <- as.numeric(outlier.mat$outliers)
  dfout <- specdf[setdiff(rownames(specdf),outlier.vector),]
        
dfout
        
}