referenz <- 35.01
datei <- read.table("foo.txt", 
               sep="\n", 
               fill=FALSE, 
               strip.white=TRUE);

datei <- matrix(t(datei), ncol=3, byrow=TRUE)
colnames(datei) <- c("id", "flags", "runtime")
temp <-gsub("[ ]+[(]sec[])]$", "", datei[,'runtime'])
temp <-as.numeric(gsub("^time[ ]+=[ ]+", "", temp))
datei[,'runtime'] <- temp

speedup <- referenz/temp
ord <- order(speedup)
speedup <- sort(speedup)
datei <- cbind(datei[ord,], speedup)

write.csv(datei, file = "foo.csv", quote = FALSE, na = "NA", row.names = FALSE)
