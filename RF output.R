testdf <- res.flow$confusion

tab <- testdf$table
#str(tab)
d <- rowSums(tab[1:nrow(tab),])#totals
c <- colSums(tab[,1:ncol(tab)])#totals

df <- rbind(tab, c)
df <- cbind(df, d)

#calculate diagonale
df[nrow(df),ncol(df)] <- 0
df[nrow(df),ncol(df)] <- sum(diag(df))

di <- diag(df)

PA <- round(head(di, -1)/head(df['c',],-1),3)*100#PA
UA <- round(head(di, -1)/head(df[,'d'],-1),3)*100#UA
