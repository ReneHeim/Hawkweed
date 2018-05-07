rmv.funct.outlier <- function(data){
        
        require(fda)
        require(fda.usc)
        
        #data <- subsets[[1]] 
        
        i <- seq(2, ncol(data))
        data.wo.noise.mat <- as.matrix(data[,i]) #As matrix to be able to transform the object (1452-52=1400)
        
        labnames <- list(main="Width", xlab="Wavelength [nm]", ylab="Reflectance [%]") #labnames to have plot information within fdata object
        myfdata <- fdata(data.wo.noise.mat, argvals = as.integer(names(data[,i])), names = labnames) #Why as integer??
        
        outlier.mat <- outliers.depth.trim(myfdata, dfunc = depth.mode, nb = 10, smo = 0.2, trim = 0.1, ns = 0.5) #Smoothing variables here are a guess
        
        outlier.vector <- as.numeric(outlier.mat$outliers)
        
        as.data.frame(data[-outlier.vector, ])
        
}