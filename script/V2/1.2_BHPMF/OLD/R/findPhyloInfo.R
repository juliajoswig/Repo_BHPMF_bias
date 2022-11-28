# hierarchy.info <- read.table(hierarchy_file_name, sep = "\t")
PreprocessHierarchy <- function(X, hierarchy.info, tmp.dir) {
num.levels <- ncol(hierarchy.info)
num.rows <- nrow(hierarchy.info)

level.names <- vector("list", num.levels)
num.nodes.per.level <- rep(0, num.levels)
for(level in 1 : num.levels) {
	level.names[[num.levels - level + 1]] <- sort(unique(hierarchy.info[, level]))
	num.nodes.per.level[num.levels - level + 1] <- length(level.names[[num.levels - level + 1]])
}
sum.num.nodes.per.level <- sum(num.nodes.per.level)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#change the string names to index
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
hierarchy.info.idx <- matrix(numeric(0), num.rows, num.levels - 1);

hierarchy.info.data.frame <- as.data.frame(hierarchy.info)
for(level in seq(1,num.levels-1, by=1)) {
	  tmp.str <- c(hierarchy.info.data.frame[, level+1])
	  hierarchy.info.idx[, level] <- c(tmp.str);
}

file.name <- paste(tmp.dir, '/processed_hierarchy_info.Rda', sep="");
save(file = file.name, list = c('hierarchy.info.idx', 'level.names', 'hierarchy.info' ,'num.nodes.per.level'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#find phylo info at higher levels
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index <- 1:num.rows
old.child.info <- c(0)

num.children <- rep(0, sum.num.nodes.per.level)
num.parent <- rep(1, sum.num.nodes.per.level)
parent <- rep(0, sum.num.nodes.per.level)

start.ind <- sum(num.nodes.per.level[1 : num.levels - 1]) + 1
last.ind <- start.ind + num.rows - 1
parent[start.ind:last.ind] <- hierarchy.info.idx[, 1]

for (level in seq(num.levels - 1, 2, by = -1)) {
	start.ind = sum(num.nodes.per.level[1:level-1])
	len = length(level.names[[level]])
	new.child.info <- rep(0, length(level.names[[level-1]]))
	# cat("level: ", level, "\n")
	for (ind in 1 : len) {
        idx <- index[hierarchy.info.idx[, num.levels-level] == ind];  # index of all children belong to this node.
		par.id <- unique(hierarchy.info.idx[idx, num.levels - level + 1])  # check if genus information is unique.
		if (length(par.id) > 1) {
			print(paste("Error in Hierarchy: node '",level.names[[leve]][ind], "' at level ", level, " belongs to more than one parent (group) :"))
			print(paste(level.names[[level-1]][par.id]))
			break
		}
	
		start.ind <- start.ind + 1
		#add number of children
		if (level == num.levels-1) {
			num.children[start.ind] <- length(idx)
		} else {
			num.children[start.ind] <- old.child.info[ind]
		}
		# cat(ind, "\t")
		# cat("par.id: " , par.id, "\n")
		num.parent[start.ind] <- 1		
		parent[start.ind] <- par.id		
		new.child.info[par.id] <- new.child.info[par.id] + 1
	}
	rm(old.child.info)
	old.child.info <- new.child.info
}

num.parent[1:num.nodes.per.level[1]] <- 0
num.children[1:num.nodes.per.level[1]] <- old.child.info

file.name = paste(tmp.dir, "/num_parents.txt", sep = "")
write.table(num.parent, file = file.name, sep = "\t", col.names = F, row.names = F)

file.name = paste(tmp.dir, "/num_children.txt", sep="")
write.table(num.children, file = file.name, sep = "\t", col.names = F, row.names = F)

file.name = paste(tmp.dir, "/hierarchy_info.txt", sep="")
write.table(parent, file = file.name, sep = "\t", col.names = F, row.names = F)
}