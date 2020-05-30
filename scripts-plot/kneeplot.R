
a <- read.table('cell_readcounts.txt', header = F, stringsAsFactors = F)
x <- cumsum(a$V1)
x <- x/max(x)
plot(1:25000, x[1:25000], type = 'l', col = 'blue')

