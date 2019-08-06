# VISTA enhancer analysis, pooled

# Load individual matrices
vista_val_pooled = lapply(as.vector(orthogonal_sets$File),function(x) read.table(paste("orthogonal/vista/",x,"_vista.txt",sep=""),sep="\t",
                                                                          col.names=c(TE_coordinates[c(1:4,6,7,5)],"chr_vista","start_vista","stop_vista","element","Validation","Overlap")))
names(vista_val_pooled) = as.vector(orthogonal_sets$File)

# Remove those without any results
vista_val_pooled = vista_val_pooled[sapply(vista_val_pooled, function(x) dim(x)[1]) > 0]

# Combine matrices
vista_val_pooled = ldply(vista_val_pooled,.id = "File")

# Add state/set information
vista_val_pooled = merge(vista_val_pooled,orthogonal_sets,by="File",all.x=TRUE)
vista_val_pooled$Set = factor(vista_val_pooled$Set,levels=names(orthogonal_labels))

# Unique VISTA enhancers overlapping TEs in that set, combining all states
vista_val_pooled_fraction = ddply(vista_val_pooled,.(Set),summarise,Total=length(unique(element)),
                           Pos=length(unique(element[which(Validation == "positive")])),
                           Neg=length(unique(element[which(Validation == "negative")])),
                           Fraction=Pos/Total)