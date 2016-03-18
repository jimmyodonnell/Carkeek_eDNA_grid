primer_summary <- read.table("/Users/jimmy.odonnell/Kelly_Lab_Big/Illumina_Data_Raw/20150717/libraries/all/primer_summary_re.txt", header = TRUE, stringsAsFactors = FALSE)

plot(primer_summary$X12SF/primer_summary$total_seq)

lib_PhiX <- grep("Undetermined", primer_summary$filename)

primer_proportion <- primer_summary[-lib_PhiX,3:6]/primer_summary[-lib_PhiX,2]

stripchart(primer_proportion, las = 1, pch = 1, method = "jitter")

contam <- 

primer_proportion[,]

which(primer_proportion[,c("X16SF", "X16SR")] < 0.5)
p12S_in_16S <- primer_proportion[,c("X12SF", "X12SR")][primer_proportion[,c("X12SF", "X12SR")] < 0.2]

range(p12S_in_16S)
round(mean(p12S_in_16S), digits = 5); round(sd(p12S_in_16S), digits = 5); round(max(p12S_in_16S), digits = 5)
primer_summary_long <- melt(primer_summary[,-2], id = c("filename"))
pdf(file = "primer_summary.pdf", width = 7, height = 3)
p <- ggplot(aes(x = filename, weight=value/1000000, fill=variable), data = primer_summary_long)
p + geom_bar() + 
  scale_fill_discrete("Primer") + 
  scale_x_discrete(breaks = primer_summary$filename, labels = 1:nrow(primer_summary)) + 
  labs(x="Libraries", y="Millions of reads", title="Primer Sequence Abundance")
dev.off()

