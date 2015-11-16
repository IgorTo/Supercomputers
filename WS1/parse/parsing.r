#referenz <- 35.01
datei <- read.table("gcc_serial_results.log", 
               sep="\n", 
               fill=FALSE, 
               strip.white=TRUE);

datei <- matrix(t(datei), ncol=3, byrow=TRUE)
colnames(datei) <- c("id", "flags", "runtime")
temp <-gsub("[ ]+[(]s[])]$", "", datei[,'runtime'])
temp <-as.numeric(gsub("^Elapsed time[ ]+=[ ]+", "", temp))
datei[,'runtime'] <- temp

nn <- nrow(datei)
referenz <- temp[nn]
datei <- datei[-nn, ]

speedup <- referenz/temp[-nn]
ord <- order(speedup)
speedup <- sort(speedup)
datei <- cbind(datei[ord,], speedup)

datei <- rbind(datei, c("0", "no flag - reference", referenz, 1))
write.csv(datei, file = "gcc.csv", quote = FALSE, na = "NA", row.names = FALSE)
