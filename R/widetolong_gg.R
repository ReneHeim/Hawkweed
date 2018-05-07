widetolong.gg <- function(datawide){
        
        datamelt<-melt(datawide, id=c("Type"))
        datamelt<-rename(datamelt, Wavelength=variable)
        datamelt<-rename(datamelt, Reflectance=value)
        
        datamelt
}