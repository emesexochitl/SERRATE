testintron <- read.table(file = "maize_intnum_drought_3h_clean_forstat.txt")
allintron <- read.table(file = "maize_intronnumv4.txt")
colnames(allintron) <- c("gene_id", "int_num")
colnames(testintron) <- c("sample", "gene_id", "int_num", "fc")
wilcox.test(allintron$int_num, testintron$int_num, correct = F)
sum(testintron$int_num)