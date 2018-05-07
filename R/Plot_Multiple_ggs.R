plot.multiple.ggs <- function(data){
        
        
        gg <- as.data.frame(data)
        ggplot(gg,aes(x=Wavelength, y=Reflectance, colour=Type))+geom_line(size=0.01)
        
}