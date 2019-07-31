# Load correlation matrices from Roadmap
cor_matrices_files = list.files(path="raw_data/correlation_matrices",pattern=".RData", full.names = TRUE)
cor_marks = lapply(cor_matrices_files,function(x) mget(load(x)))
names(cor_marks) = gsub(".RData","",gsub("raw_data/correlation_matrices/cor_","",cor_matrices_files))

# Reformat
cor_marks = lapply(cor_marks,function(x) {y = melt(as.matrix(x$markcor)); return(y)})
cor_marks = ldply(cor_marks,.id="Mark")
colnames(cor_marks)[2:4] = c("Sample1","Sample2","Correlation")

# Remove self correlations
cor_marks = cor_marks[which(cor_marks$Sample1 != cor_marks$Sample2),]

# Correlations where Sample 1 is in a replicate pair
cor_marks_replicates = cor_marks[which(cor_marks$Sample1 %in% as.vector(unique(replicates$`Sample 1`))),]

# Assign rank to each correlation
cor_marks_replicates = ddply(cor_marks_replicates,.(Sample1,Mark),transform,Rank=rank(Correlation))

# Identify replicate/different tissue pairs
cor_marks_replicates = merge(cor_marks_replicates,replicates[,c("Sample 1","Sample 2","Set")],
                             by.x=c("Sample1","Sample2"),by.y=c("Sample 1","Sample 2"),all.x=TRUE)
cor_marks_replicates[is.na(cor_marks_replicates)] = "Other"
cor_marks_replicates$Set = factor(cor_marks_replicates$Set,levels=c("Replicate","Tissue","Other"))

# Plot
ggplot(cor_marks_replicates,aes(x=Sample1,y=Correlation,color=Set)) + geom_jitter(width = 0.2) + 
  facet_wrap(~Mark) + scale_color_manual(values=c("red","blue","grey"),name="Pair type") + xlab("Sample 1")
