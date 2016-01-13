file = "../Climate_90KSNP_1000g_FinalDataMatrix.csv"

d = read.csv(file)

#sub_d = d[,c(1,2,23,24,25,12,19)]

#write.csv(sub_d, file = "extracted_phenotype.csv", row.names=F)

markers = colnames(d)[c(1, 82:ncol(d))]
sub_d = d[, c(1, 82:ncol(d))]
sub_d = t(sub_d)
rownames(sub_d) = markers

write.csv(sub_d, file="extracted_genotype.csv")
