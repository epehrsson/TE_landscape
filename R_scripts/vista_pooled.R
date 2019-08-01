# VISTA enhancer analysis, pooled

# Load individual matrices
vista_val = lapply(as.vector(orthogonal_sets$File),function(x) read.table(paste("orthogonal/vista/",x,"_vista.txt",sep=""),sep="\t",
                                                                          col.names=c(TE_coordinates[c(1:4,6,7,5)],"chr_vista","start_vista","stop_vista","element","Validation","Overlap")))
names(vista_val) = as.vector(orthogonal_sets$File)

# Remove those without any results
vista_val = vista_val[sapply(vista_val, function(x) dim(x)[1]) > 0]

# Combine matrices
vista_val = ldply(vista_val,.id = "File")

# Add state/set information
vista_val = merge(vista_val,orthogonal_sets,by="File",all.x=TRUE)
vista_val$Set = factor(vista_val$Set,levels=names(orthogonal_labels))

# Unique VISTA enhancers overlapping TEs in that set, combining all states
vista_val_fraction = ddply(vista_val,.(Set),summarise,Total=length(unique(element)),
                           Pos=length(unique(element[which(Validation == "positive")])),
                           Neg=length(unique(element[which(Validation == "negative")])),
                           Fraction=Pos/Total)