outlier.test <- function(data.before, data.after){
        
        removed=as.numeric(nrow(data.before)-nrow(data.after))
        
        message(sprintf("Outlier removed: %s\n", removed))
}