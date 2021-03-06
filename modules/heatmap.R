#!/usr/bin/Rscript
library(graphics)
x11()
file_name = (commandArgs(TRUE)[1])
out_file = (commandArgs(TRUE)[2])
dsRNA_targ = read.csv(file_name, 
                        header = FALSE, 
                        row.names=1, colClasses=c("numeric","numeric"))
z=c(dsRNA_targ[,1])
z=matrix(z)
my_palette <- colorRampPalette(c("red","black","purple","green"))(n = 1000)
image(z,col = my_palette, axes = FALSE)
xspace=seq(0,1,(250/length(z)))
xlabels = seq(0,length(z),250)
axis(1, at= xspace, labels=xlabels,las=2)
dev.copy(jpeg, filename="plot.jpg")
dev.off()