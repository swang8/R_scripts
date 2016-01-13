
### Manhattan plot
files <- commandArgs(trailingOnly = T)
accession = gsub("(.*).GWAS.*", "\\1", files[1]);
data <- read.table(files[1], header=T, sep=",")
## SNP,Chromosome,Position ,P.value
data <- data[,c(1,1,2,3,4)]
data[,5] = -1 * log(data[,5]) / log(10)
colnames(data) <- c("Marker", "SNP", "Chr", "Position", accession)
data <- data[order(data$Chr), ]

#chrs <- paste(rep(1:7, each=6), rep(c("A", "B", "D"), each=2), sep="")
#data[,3] = chrs[data[,3]]

chrs = sort(data[, 3])

genomes=c("A", "B", "D")


file = paste(accession, "man__.pdf", sep="_")
pdf(file)  # generate PDF files
par(mfrow=c(3,1))


y_max = max(data[,5]);
if (y_max < 8) y_max = 8;


for (gn in genomes)
{
  gn_chrs <- unique(chrs[grepl(gn, chrs)])
  subdata <- subset(data, grepl(gn, Chr) & (!grepl("NA", Chr)))
  for (i in gn_chrs){
   ss <- subset(subdata, grepl(i, Chr));
   ss <- ss[order(ss$Position), ]
   ss$Position = 1:length(ss$Position)
   subdata[grepl(i, subdata$Chr), ] = ss 
  }
  # test
    if(grepl("A", gn)){
      print(subset(subdata, grepl("wsnp_BM140362A_Ta_2_2", SNP)))
    }
  #
  pp <- subdata$Position
  f <- function(x){if(x<0){return(0)}else{return(x)}}
  subdata$Position <- sapply(pp,f)
  chr_col = character(length(subdata$Chr))
  for (i in 1:length(chr_col))
  {
    chr_num = as.numeric(substr(subdata$Chr[i],1,1))
    if(chr_num%%2 != 0){chr_col[i] = "slateblue"}
    if(chr_num%%2 == 0){chr_col[i] = "brown1"}
  }
 
  max_vals=numeric()
  min_vals = numeric()
  for (c in gn_chrs){
    mm <- max(subdata$Position[grepl(c, subdata$Chr)])
    mn <- min(subdata$Position[grepl(c, subdata$Chr)])
    max_vals = c(max_vals, mm)
    min_vals = c(min_vals, mn)
  }
  addit_vals <- numeric(7);
  addit_vals[1] = -1*min_vals[1];
  for (nn in 2:7)
  {
    addit_vals[nn] = sum(max_vals[1:(nn-1)]) - sum(min_vals[1:nn]) 
  }
  
  new_pos = numeric()
  ft = 1;
  for (c in 1:length(gn_chrs))
  {
    p <- subdata$Position[grepl(gn_chrs[c], subdata$Chr)]
    p <- (p+addit_vals[c])*ft
    new_pos <- c(new_pos, p)
  }
  
  l <- colnames(subdata);
  y_uplim <- max(subdata[,5:length(l)])
  if(y_uplim < 8){y_uplim = 8}
  lab_pos = (min_vals + max_vals + 2*addit_vals)*ft/2
  
  ##plot(new_pos, subdata[, 5], col=chr_col, axes=F, ylim=c(0,y_max), pch=20, xlab="", ylab="-Log10(P)", main=paste(gn, "genome"))
  plot(new_pos, subdata[, 5], col=chr_col, axes=F, ylim=c(0,y_uplim), pch=20, xlab="", ylab="-Log10(P)", main=paste(gn, "genome"))
  abline(h=5, lty=3, col="grey")
  axis(1, at=lab_pos, labels=gn_chrs, tick=F)
  axis(2)
  box()
    
}
dev.off()

