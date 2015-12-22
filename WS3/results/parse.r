hasw <- as.matrix(read.table("3_haswell.txt", 
                    skip = 8,
                    sep="\n", 
                    fill=FALSE, 
                    strip.white=TRUE))

sandy <- as.matrix(read.table("3_sandybridge.txt",
                    skip = 8,
                    sep="\n", 
                    fill=FALSE, 
                    strip.white=TRUE))
hasw <- hasw[!grepl("^Tue Dec 22 [0-9]+:[0-9]+:[0-9]+ CET 2015", hasw)]
sandy <- sandy[!grepl("^Tue Dec 22 [0-9]+:[0-9]+:[0-9]+ CET 2015", sandy)]

hasw <- matrix(hasw, ncol=3, byrow=TRUE)
sandy <- matrix(sandy, ncol=3, byrow=TRUE)
colnames(hasw) <- c("size", "computation", "mpi")
colnames(sandy) <- c("size", "computation", "mpi")


tempCh <-as.numeric(gsub("^Computation[ ]+time:[ ]+", "", hasw[,'computation']))
tempCs <-as.numeric(gsub("^Computation[ ]+time:[ ]+", "", sandy[,'computation']))
tempMh <-as.numeric(gsub("^MPI[]+time:[ ]+", "", hasw[,'mpi']))
tempMs <-as.numeric(gsub("^MPI[]+time:[ ]+", "", sandy[,'mpi']))

hasw[,'computation'] <- tempCh
hasw[,'mpi'] <- tempMh
sandy[,'computation'] <- tempCs
sandy[,'mpi'] <- tempMs

rowS <- 1:nrow(sandy)
rowH <- 1:nrow(hasw)

avgSandyC <- c(1,1,1,1,1,1,1)
avgHaswC <- c(1,1,1,1,1,1,1)
avgSandyM <- c(1,1,1,1,1,1,1)
avgHaswM <- c(1,1,1,1,1,1,1)

for (i in 1:6){
  avgSandyC[i] <- sum(as.numeric(sandy[ rowS%%7 == i, "computation"]))/10
  avgSandyM[i] <- sum(as.numeric(sandy[ rowS%%7 == i, "mpi"]))/10
  avgHaswC[i] <- sum(as.numeric(hasw[ rowH%%7 == i, "computation"]))/10
  avgHaswM[i] <- sum(as.numeric(hasw[ rowH%%7 == i, "mpi"]))/10
}
avgSandyC[7] <- sum(as.numeric(sandy[ rowS%%7 == 0, "computation"]))/10
avgSandyM[7] <- sum(as.numeric(sandy[ rowS%%7 == 0, "mpi"]))/10
avgHaswC[7] <- sum(as.numeric(hasw[ rowH%%7 == 0, "computation"]))/10
avgHaswM[7] <- sum(as.numeric(hasw[ rowH%%7 == 0, "mpi"]))/10

sizes <- c("64", "128", "256", "512", "1024", "2048", "4069")



#-------------------------SANDY

datei1 <- cbind(sizes, avgSandyC, avgSandyM)
#slikca

colnames(datei1) <- c("size", "avg.computation", "avg.mpi")
write.csv(datei1, file = "sandy.csv", quote = FALSE, na = "NA", row.names = FALSE)


#-------------------------HASWELL

datei2 <- cbind(sizes, avgHaswC, avgHaswM) 
#slikca

colnames(datei2) <- c("size", "avg.computation", "avg.mpi")
write.csv(datei2, file = "haswell.csv", quote = FALSE, na = "NA", row.names = FALSE)




