testdf <- res.flow$confusion

tab <- testdf$table
#str(tab)
d <- rowSums(tab[1:10,])
c <- colSums(tab[,1:10])

df <- rbind(tab, c)
df <- cbind(df, d)

#calculate diagonale
df[11,11] <- 0
df[11,11] <- sum(diag(df))
