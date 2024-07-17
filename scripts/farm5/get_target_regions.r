library(stringr)
library(stringi)
library(Biostrings)
library(parallel)

arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

threshold.variance <- 0.2

sample.list <- arg.values[1]
sample.names <- setNames(nm = read.table(sample.list, stringsAsFactors = F)[[1]])

gene.regions.file <- arg.values[2]
all.gene.coordinates <- read.table(gene.regions.file, sep = '\t', header = T, row.names = 1)

meta.file <- arg.values[3]
meta <- read.table(meta.file, sep = '\t', header = T, row.names = 1, quote = '', comment.char = '')[sample.names, ]
expected.copy.number.on.X <- c(2, 1)[(meta$sex_call == 'M') + 1]

cov.var.file <- arg.values[4]
cov.var <- read.table(cov.var.file, header = T, sep = '\t', row.names = 1)
high.var.samples <- intersect(sample.names, rownames(cov.var)[cov.var$autosomes > threshold.variance])

coverage.folder <- arg.values[5]
diagnostic.reads.folder <- arg.values[6]
plotting.functions.file <- arg.values[7]
num.cores <- as.numeric(arg.values[8])

# Function to load up each file and store the results in a list
load.hmm.file <- function(sample.name, folder = '.'){
	all.files <- list.files(folder)
	# Find the file corresponding to this sample
	if (length(grep(sample.name, all.files, value = T)) == 0)
		stop(paste('No matching files were found for sample', sample.name))
	this.file <- paste(folder, grep(sample.name, all.files, value = T), sep = '/')
	if (length(this.file) > 1)
		stop(paste('More than one matching file were found for sample', sample.name))
	# Load the file 
	cat('\t', this.file, '\n', sep = '')
	read.table(this.file, header = T)
}

# Function that will take a table of HMM output and shrink it to contain only a desired region
shrink.compact.hmm <- function(region.coords, input.hmm){
	subset(input.hmm, (Position >= region.coords[1]) & (Position <= region.coords[2]))
}

# In order to process each hmm file only once, our function for multi-threading needs to output all 
# objects needed from each hmm file
process_2RL_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, '2RL/HMM_output', sep = '/'))
	output.list <- list(Cyp6 = shrink.compact.hmm(region.coords$Cyp6, hmm.table),
						Ache = shrink.compact.hmm(region.coords$Ache, hmm.table),
						Esterase = shrink.compact.hmm(region.coords$Esterase, hmm.table),
						Gst1 = shrink.compact.hmm(region.coords$Gst1, hmm.table),
						Harb1 = shrink.compact.hmm(region.coords$Harb1, hmm.table),
						Gst2 = shrink.compact.hmm(region.coords$Gst2, hmm.table),
						Zinccarbo = shrink.compact.hmm(region.coords$Zinccarbo, hmm.table))
	return(output.list)
}

process_3RL_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, '3RL/HMM_output', sep = '/'))
	output.list <- list(GABA = shrink.compact.hmm(region.coords$GABA, hmm.table),
						Carboxypep = shrink.compact.hmm(region.coords$Carboxypep, hmm.table))
	return(output.list)
}

process_X_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, 'X/HMM_output', sep = '/'))
	output.list <- list(Cyp9 = shrink.compact.hmm(region.coords$Cyp9, hmm.table),
						RDGA = shrink.compact.hmm(region.coords$RDGA, hmm.table),
						NADHCyp = shrink.compact.hmm(region.coords$NADHCyp, hmm.table))
	return(output.list)
}


# Function that will calculate the mode of the HMM coverage in a given region
get.gene.mode <- function(hmm.data, target.region, window.size = 300){
	target.region <- as.numeric(target.region)
	hmm.in.region <- subset(hmm.data, Position >= (target.region[1] - window.size) & Position < target.region[2])$CNV
	# We calculate the mode. Where there is more than one mode, we take the smallest one
	un <- unique(hmm.in.region)
	tab <- tabulate(match(hmm.in.region, un))
	hmm.mode <- min(un[tab == max(tab)])
	hmm.mode
}

# For funestus, replace this with a list of genes inside the region that you want to display (for Cyp6, 
# I believe we did 8300000:9200000). Here I have just picked a couple of genes at random and given them
# random names. 
focal.genes <- list(Cyp6 = c(    nadh =         'LOC125764767',
                                 Cyp6a1 =       'LOC125764719',
                                 Cyp6a2 =       'LOC125764726',
                                 Cyp6a3 =       'LOC125764717',
                                 Cyp6a4 =       'LOC125764715',
                                 Cyp6a5 =       'LOC125764705',
                                 Chan1  =       'LOC125761428',
                                 Chan2  =       'LOC125764677',
                                 estern =       'LOC125764700',
                                 Cyp6n1 =       'LOC125764711',
                                 Cyp6n2 =       'LOC125764713',
                                 Chan3  =       'LOC125764682',
                                 Cyp6a6 =       'LOC125764714',
                                 Cyp6a7 =       'LOC125764712'),
                    Ache = c(    Ache =         'LOC125763056'),					
                    Esterase = c(AGAP003155 =   'LOC125762085',
                                 aldehyd =      'LOC125762067'),
                    Gst1 = c(    Gstd1 =        'LOC125761308',
                                 Gstd2 =        'LOC125761316',
                                 Gstd3 =        'LOC125761314',
                                 Gstd4 =        'LOC125761303',
                                 Gstd5 =        'LOC125761313',
                                 Gstd6 =        'LOC125761315'),
	                Harb1 = c(   Harb12 =       'LOC125762969',
                                 Harb11 =       'LOC125760979'),
                    Gst2 = c(    gste1 =        'LOC125763984',
                                 gste2 =        'LOC125763977',
                                 gste3 =        'LOC125763983',
                                 gste4 =        'LOC125763985',
                                 gste5 =        'LOC125763986',
                                 gste6 =        'LOC125763982'),
		            Zinccarbo =c(Zinccarbo1 =   'LOC125763331',
                                 Zinccarbo2 =   'LOC125763332',
                                 Zinccarbo3 =   'LOC125763330',
                                 Zinccarbo4 =   'LOC125763333'),
                    GABA = c(    Mono5 =        'LOC125770368',
                                 Mono3 =        'LOC125769985',
                                 GABA =         'LOC125769835'),
                    Carboxypep=c(Carboxy1 =     'LOC125767895',
                                 Carboxy2 =     'LOC125767890',
                                 Carboxy3 =     'LOC125767907',
                                 Carboxy4 =     'LOC125767894',
                                 Carboxy5 =     'LOC125767904',
                                 Carboxy6 =     'LOC125767898'),
                    Cyp9 = c(    cuticle =      'LOC125764385',
                                 Cyp1 =         'LOC125764232'),
		            RDGA = c(    rdga =	        'LOC125769763',
                                 abc  =	        'LOC125760949'),
                    NADHCyp = c( nadcyp =       'LOC125760701',
                                 acbe =         'LOC125769763')
)

# Get the genetic coordinates of those genes.
gene.coords <- lapply(focal.genes, function(x) {X = all.gene.coordinates[x, ]; rownames(X) = names(x); X})

# Define the coordinates of the genetic regions that we will use for each cluster
region.coords = list(Cyp6 = c(8240000, 9820000),
                     Ache = c(19400000, 19500000), 
                     Esterase = c(39020000, 39850000), 
                     Gst1 = c(42000000, 42460000), 
                     Harb1 = c(57054400, 58000000), 
                     Gst2 = c(76000000, 76900000),
                     Zinccarbo = c(85450000, 85500000), 
                     GABA = c(13400000, 14200000), 
                     Carboxypep = c(30300000, 30350000), 
                     Cyp9 = c(8350000, 8870000), 
                     RDGA = c(13590000, 14000000),
                     NADHCyp = c(14350000, 14600000)
					)

plotting.ranges <- list()
compact.hmm <- list()
cov.cnv.samples <- list()
# Create an empty dataframe where each row is a sample
hmm.cnv.table <- data.frame(row.names = sample.names)
# We will include a table giving the sex of the sample, since this is relevant to copy number of any region on the X
hmm.cnv.table$sex <- meta[rownames(hmm.cnv.table), 'sex_call']
# And a column showing whether the sample has high coverage variance, since these will be dubious
hmm.cnv.table$high.var <- rownames(hmm.cnv.table) %in% high.var.samples


#### 2RL
# Load up the 2RL data and get the tables for the Ace1 and Cyp6 regions
cat('Creating shrunk table for regions on 2RL.\n')
temp.2RL <- mclapply(sample.names, process_2RL_regions, mc.cores = num.cores)

# Because the above line uses parallel processing, I am not 100% confident that it will correctly preserve 
# sample order, so in the next couple of lines we use lapply on the sample.names vector to make sure the order 
# is correct. 
compact.hmm$Cyp6 <- lapply(sample.names, function(s) temp.2RL[[s]]$Cyp6)
compact.hmm$Ache <- lapply(sample.names, function(s) temp.2RL[[s]]$Ache)
compact.hmm$Esterase <- lapply(sample.names, function(s) temp.2RL[[s]]$Esterase)
compact.hmm$Gst1 <- lapply(sample.names, function(s) temp.2RL[[s]]$Gst1)
compact.hmm$Harb1 <- lapply(sample.names, function(s) temp.2RL[[s]]$Harb1)
compact.hmm$Gst2 <- lapply(sample.names, function(s) temp.2RL[[s]]$Gst2)
compact.hmm$Zinccarbo <- lapply(sample.names, function(s) temp.2RL[[s]]$Zinccarbo)

rm(temp.2RL) 
# Get the ranges within this region that we wish to use for plotting (here I've just made it the same as the 
# total region in which we looked for discordant reads
plotting.ranges$Cyp6 <- region.coords$Cyp6
plotting.ranges$Ache <- region.coords$Ache
plotting.ranges$Esterase <- region.coords$Esterase
plotting.ranges$Gst1 <- region.coords$Gst1
plotting.ranges$Harb1 <- region.coords$Harb1
plotting.ranges$Gst2 <- region.coords$Gst2
plotting.ranges$Zinccarbo <- region.coords$Zinccarbo

# Record the CNVs that were discovered based on coverage in these regions. Here we list the genes that we want
# to have modal-coverage data for
# Ace1
hmm.cnv.table[['Ache']] <- unlist(lapply(compact.hmm$Ache, get.gene.mode, target.region = gene.coords$Ache['Ache', 2:3])) - 2
cov.cnv.samples$Ache <- rownames(subset(hmm.cnv.table, Ache > 0))

#Cyp6
cyp6.genes <- c('Cyp6a1', 'Cyp6a2', 'Cyp6a3', 'Cyp6a4', 'Cyp6a5', 'Cyp6a6', 'Cyp6a7', 'Cyp6n1', 'Cyp6n1')
for (g in cyp6.genes)
	hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Cyp6, get.gene.mode, target.region = gene.coords$Cyp6[g, 2:3])) - 2
hmm.cnv.table[["Max_Cyp6"]] <- apply(hmm.cnv.table[, cyp6.genes], 1, max)
cov.cnv.samples$Cyp6 <- rownames(subset(hmm.cnv.table, Max_Cyp6 > 0))

#Esterase
for (g in names(focal.genes$Esterase))
	hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Esterase, get.gene.mode, target.region = gene.coords$Esterase[g, 2:3])) - 2
hmm.cnv.table[["Max_Esterase"]] <- apply(hmm.cnv.table[, names(focal.genes$Esterase)], 1, max)
cov.cnv.samples$Esterase <- rownames(subset(hmm.cnv.table, Max_Esterase > 0))

#Gst1
for (g in names(focal.genes$Gst1))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Gst1, get.gene.mode, target.region = gene.coords$Gst1[g, 2:3])) - 2
hmm.cnv.table[["Max_Gst1"]] <- apply(hmm.cnv.table[, names(focal.genes$Gst1)], 1, max)
cov.cnv.samples$Gst1 <- rownames(subset(hmm.cnv.table, Max_Gst1 > 0))

#Harb1
for (g in names(focal.genes$Harb1))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Harb1, get.gene.mode, target.region = gene.coords$Harb1[g, 2:3])) - 2
hmm.cnv.table[["Max_Harb1"]] <- apply(hmm.cnv.table[, names(focal.genes$Harb1)], 1, max)
cov.cnv.samples$Harb1 <- rownames(subset(hmm.cnv.table, Max_Harb1 > 0))

#Gst2
for (g in names(focal.genes$Gst2))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Gst2, get.gene.mode, target.region = gene.coords$Gst2[g, 2:3])) - 2
hmm.cnv.table[["Max_Gst2"]] <- apply(hmm.cnv.table[, names(focal.genes$Gst2)], 1, max)
cov.cnv.samples$Gst2 <- rownames(subset(hmm.cnv.table, Max_Gst2 > 0))

#Zinccarbo
for (g in names(focal.genes$Zinccarbo))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Zinccarbo, get.gene.mode, target.region = gene.coords$Zinccarbo[g, 2:3])) - 2
hmm.cnv.table[["Max_Zinccarbo"]] <- apply(hmm.cnv.table[, names(focal.genes$Zinccarbo)], 1, max)
cov.cnv.samples$Zinccarbo <- rownames(subset(hmm.cnv.table, Max_Zinccarbo > 0))


#### 3RL
# Load up the 3RL data and get the tables for the GABA(rdl) and Carboxyprep regions
cat('Creating shrunk table for regions on 3RL.\n')
temp.3RL <- mclapply(sample.names, process_3RL_regions, mc.cores = num.cores)

compact.hmm$GABA <- lapply(sample.names, function(s) temp.3RL[[s]]$GABA)
compact.hmm$Carboxypep <- lapply(sample.names, function(s) temp.3RL[[s]]$Carboxypep)
rm(temp.3RL)

# Get the ranges within this region that we wish to use for plotting (here I've just made it the same as the 
# total region in which we looked for discordant reads
plotting.ranges$GABA <- region.coords$GABA
plotting.ranges$Carboxypep <- region.coords$Carboxypep

# Record the CNVs that were discovered based on coverage in these regions. Here we list the genes that we want
# to have modal-coverage data for
#GABA
for (g in names(focal.genes$GABA))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$GABA, get.gene.mode, target.region = gene.coords$GABA[g, 2:3])) - 2
hmm.cnv.table[["Max_GABA"]] <- apply(hmm.cnv.table[, names(focal.genes$GABA)], 1, max)
cov.cnv.samples$GABA <- rownames(subset(hmm.cnv.table, Max_GABA > 0))

#Carboxypep
for (g in names(focal.genes$Carboxypep))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Carboxypep, get.gene.mode, target.region = gene.coords$Carboxypep[g, 2:3])) - 2
hmm.cnv.table[["Max_Carboxypep"]] <- apply(hmm.cnv.table[, names(focal.genes$Carboxypep)], 1, max)
cov.cnv.samples$Carboxypep <- rownames(subset(hmm.cnv.table, Max_Carboxypep > 0))

#### X
# Load up the X data and get the tables for the Cyp9 and NADHCyp regions
cat('Creating shrunk table for regions on X.\n')
temp.X <- mclapply(sample.names, process_X_regions, mc.cores = num.cores)

# Because the above line uses parallel processing, I am not 100% confident that it will correctly preserve 
# sample order, so in the next couple of lines we use lapply on the sample.names vector to make sure the order 
# is correct.
compact.hmm$Cyp9 <- lapply(sample.names, function(s) temp.X[[s]]$Cyp9)
compact.hmm$RDGA <- lapply(sample.names, function(s) temp.X[[s]]$RDGA)
compact.hmm$NADHCyp <- lapply(sample.names, function(s) temp.X[[s]]$NADHCyp)
rm(temp.X)

# Get the ranges within this region that we wish to use for plotting (here I've just made it the same as the 
# total region in which we looked for discordant reads
plotting.ranges$Cyp9 <- region.coords$Cyp9
plotting.ranges$RDGA <- region.coords$RDGA
plotting.ranges$NADHCyp <- region.coords$NADHCyp

# Record the CNVs that were discovered based on coverage in these regions. Here we list the genes that we want
# to have modal-coverage data for
#Cyp9
for (g in names(focal.genes$Cyp9))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Cyp9, get.gene.mode, target.region = gene.coords$Cyp9[g, 2:3])) - 2
hmm.cnv.table[["Max_Cyp9"]] <- apply(hmm.cnv.table[, names(focal.genes$Cyp9)], 1, max)
cov.cnv.samples$Cyp9 <- rownames(subset(hmm.cnv.table, Max_Cyp9 > 0))

#RDGA
for (g in names(focal.genes$RDGA))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$RDGA, get.gene.mode, target.region = gene.coords$RDGA[g, 2:3])) - 2
hmm.cnv.table[["Max_RDGA"]] <- apply(hmm.cnv.table[, names(focal.genes$RDGA)], 1, max)
cov.cnv.samples$RDGA <- rownames(subset(hmm.cnv.table, Max_RDGA > 0))

#NADHCyp
for (g in names(focal.genes$NADHCyp))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$NADHCyp, get.gene.mode, target.region = gene.coords$NADHCyp[g, 2:3])) - 2
hmm.cnv.table[["Max_NADHCyp"]] <- apply(hmm.cnv.table[, names(focal.genes$NADHCyp)], 1, max)
cov.cnv.samples$NADHCyp <- rownames(subset(hmm.cnv.table, Max_NADHCyp > 0))


# Get a named vector of the types of discordant reads that will be used 
discordant.read.types <- c(FA = 'FA', SS = 'SS', FM = 'FM', XC = 'crosschrom')

# Write a function to load the results of discordant read analysis from file for a given sample
get.discordant.reads <- function(this.sample, target.folder, target.region.coords){
	this.file <- grep(paste(this.sample, '.*csv$', sep = ''), list.files(target.folder), value = T)
	if (length(this.file) != 1)
		stop('There should be one and only one file in the folder that matches the sample name.')
   	pos.table <- read.table(paste(target.folder, this.file, sep='/'), header = T, sep='\t', colClasses = c('character', 'integer', 'character'))
	split.pos.table <- split(pos.table, pos.table$Type)
	get.discordant.read.type <- function(read.type){
		full.read.type <- paste(read.type, 'mapq >= 10')
		if (!(full.read.type %in% names(split.pos.table))){
			return(data.frame(Position = integer(), Mate.position = integer()))
		}
		allpos <- split.pos.table[[full.read.type]][, c('Position', 'Mate.position')]
		if (read.type %in% c('FA', 'SS', 'FM')){
			allpos$Mate.position <- as.integer(allpos$Mate.position)
			which.in.region <- ((allpos$Position >= target.region.coords[1]) & (allpos$Position <= target.region.coords[2])) | ((allpos$Mate.position >= target.region.coords[1]) & (allpos$Mate.position <= target.region.coords[2]))
		}
		else{
			which.in.region <- (allpos$Position >= target.region.coords[1]) & (allpos$Position <= target.region.coords[2])
		}
		allpos <- allpos[which.in.region, ]
		allpos[order(allpos$Position),]
	}
	lapply(discordant.read.types, get.discordant.read.type)
}

# Write a function to load the results of breakpoint read analysis from file for a given sample
get.breakpoint.reads <- function(this.sample, target.folder, target.region.coords){
	this.file <- grep(paste(this.sample, '.*csv$', sep = ''), list.files(target.folder), value = T)
	if (length(this.file) != 1)
		stop('There should be one and only one file in the folder that matches the sample name.')
   	pos.table <- read.table(paste(target.folder, this.file, sep='/'), header = T, sep='\t', colClasses = c('character', 'integer', 'character'))
	which.in.region <- (pos.table$Position >= target.region.coords[1]) & (pos.table$Position <= target.region.coords[2])
	pos.table <- pos.table[which.in.region, ]
	clipping.start.point <- pos.table[pos.table$Type == 'soft clipping start point mapq >= 10', c('Position', 'Clipped_sequence')]
	clipping.end.point <- pos.table[pos.table$Type == 'soft clipping end point mapq >= 10', c('Position', 'Clipped_sequence')]
	#
	list(CSP = clipping.start.point[order(clipping.start.point$Position), ],
	     CEP = clipping.end.point[order(clipping.end.point$Position), ])
}

get.diagnostic.reads <- function(this.sample, disc.target.folder, bp.target.folder, target.region.coords){
	cat('\tLoading diagnostic reads for ', this.sample, '.\n', sep = '')
	output <- get.discordant.reads(this.sample, disc.target.folder, target.region.coords)
	output[['BP']] <- get.breakpoint.reads (this.sample, bp.target.folder, target.region.coords)
	output
}

# Write a function to load the results of discordant and breakpoint reads from file for a vector of samples
get.diagnostic.reads.allsamples <- function(sample.names, disc.target.folder, bp.target.folder, target.region.coords){
	cat('Loading discordant reads from ', disc.target.folder, ' and breakpoint reads from ', bp.target.folder, '.\n', sep = '')
	lapply(sample.names, get.diagnostic.reads, disc.target.folder, bp.target.folder, target.region.coords)
}

# Set the folders in which to look for FA, SS, FM and XC reads. 
SSFA.folders <- setNames(paste(diagnostic.reads.folder, c('SSFA/2RL/Cyp6_region',
							  'SSFA/2RL/Ache_region',
							  'SSFA/2RL/Esterase_region',
							  'SSFA/2RL/Gst1_region',
							  'SSFA/2RL/Harb1_region',
							  'SSFA/2RL/Gst2_region',
							  'SSFA/2RL/Zinccarbo_region',
							  'SSFA/3RL/GABA_region',
							  'SSFA/3RL/Carboxypep_region',
							  'SSFA/X/Cyp9_region',
							  'SSFA/X/RDGA_region',
							  'SSFA/X/NADHCyp_region'),
                                                        sep = '/'),
                         c('Cyp6', 'Ache', 'Esterase', 'Gst1', 'Harb1', 'Gst2','Zinccarbo', 
						 'GABA', 'Carboxypep', 'Cyp9', 'RDGA', 'NADHCyp')
                        )

# Set the folders in which to look for breakpoint reads. 
breakpoints.folders <- setNames(paste(diagnostic.reads.folder, c('breakpoints/2RL/Cyp6_region',
								 'breakpoints/2RL/Ache_region',
								 'breakpoints/2RL/Esterase_region',
								 'breakpoints/2RL/Gst1_region',
								 'breakpoints/2RL/Harb1_region',
								 'breakpoints/2RL/Gst2_region',
								 'breakpoints/2RL/Zinccarbo_region',
								 'breakpoints/3RL/GABA_region',
								 'breakpoints/3RL/Carboxypep_region',
								 'breakpoints/X/Cyp9_region',
								 'breakpoints/X/RDGA_region',	
								 'breakpoints/X/NADHCyp_region'),
                                                               sep = '/'),
                                c('Cyp6', 'Ache', 'Esterase', 'Gst1', 'Harb1', 'Gst2','Zinccarbo', 
								'GABA', 'Carboxypep', 'Cyp9', 'RDGA', 'NADHCyp')
                               )

cat('Loading discordant and breakpoint reads\n')
diagnostic.reads <- mcmapply(get.diagnostic.reads.allsamples,
							 SSFA.folders, breakpoints.folders, 
							 region.coords, MoreArgs = list(sample.names = sample.names),
							SIMPLIFY = F, mc.preschedule = F, mc.cores = num.cores)

save.image('target_regions_analysis/target_regions.Rdata')
source(plotting.functions.file)
save.image('target_regions_analysis/target_regions.Rdata')

