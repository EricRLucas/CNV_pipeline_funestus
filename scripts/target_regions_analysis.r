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
meta <- read.table(meta.file, sep = '\t', header = T, row.names = 1, quote = '"', comment.char = '', colClasses = c(sex_call = 'factor'))[sample.names, ]
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

focal.genes <- list(Ache = c(    Ache =         'LOC125763056'),					
					Cyp6 = c(    nadh =         'LOC125764767',
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
                    Esterase = c(AGAP003155 =   'LOC125762085',
                                 aldehyd =      'LOC125762067'),
                    Gst1 = c(    Gstd1 =        'LOC125761308',
                                 Gstd2 =        'LOC125761316',
                                 Gstd3 =        'LOC125761314',
                                 Gstd4 =        'LOC125761303',
                                 Gstd5 =        'LOC125761313',
                                 Gstd6 =        'LOC125761315'),
                    Gst2 = c(    gste1 =        'LOC125763984',
                                 gste2 =        'LOC125763977',
                                 gste3 =        'LOC125763983',
                                 gste4 =        'LOC125763985',
                                 gste5 =        'LOC125763986',
                                 gste6 =        'LOC125763982'),
	                Harb1 = c(   Harb12 =       'LOC125762969',
                                 Harb11 =       'LOC125760979'),
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
region.coords = list(Ache = c(19400000, 19500000), 
                     Cyp6 = c(8240000, 9820000),
                     Esterase = c(39020000, 39850000), 
                     Gst1 = c(42000000, 42460000), 
                     Gst2 = c(76000000, 76900000),
                     Harb1 = c(57054400, 58000000), 
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

# Load up the 2RL data
cat('Creating shrunk table for regions on 2RL.\n')
temp.2RL <- mclapply(sample.names, process_2RL_regions, mc.cores = num.cores)

# Because the above line uses parallel processing, I am not 100% confident that it will correctly preserve 
# sample order, so in the next couple of lines we use lapply on the sample.names vector to make sure the order 
# is correct. 
compact.hmm$Ache <- lapply(sample.names, function(s) temp.2RL[[s]]$Ache)
compact.hmm$Cyp6 <- lapply(sample.names, function(s) temp.2RL[[s]]$Cyp6)
compact.hmm$Esterase <- lapply(sample.names, function(s) temp.2RL[[s]]$Esterase)
compact.hmm$Gst1 <- lapply(sample.names, function(s) temp.2RL[[s]]$Gst1)
compact.hmm$Gst2 <- lapply(sample.names, function(s) temp.2RL[[s]]$Gst2)
compact.hmm$Harb1 <- lapply(sample.names, function(s) temp.2RL[[s]]$Harb1)
compact.hmm$Zinccarbo <- lapply(sample.names, function(s) temp.2RL[[s]]$Zinccarbo)

rm(temp.2RL) 

# Get the ranges within this region that we wish to use for plotting (here I've just made it the same as the 
# total region in which we looked for discordant reads
plotting.ranges$Ache <- region.coords$Ache
plotting.ranges$Cyp6 <- region.coords$Cyp6
plotting.ranges$Esterase <- region.coords$Esterase
plotting.ranges$Gst1 <- region.coords$Gst1
plotting.ranges$Gst2 <- region.coords$Gst2
plotting.ranges$Harb1 <- region.coords$Harb1
plotting.ranges$Zinccarbo <- region.coords$Zinccarbo

# Record the CNVs that were discovered based on coverage in these regions
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

#Gst2
for (g in names(focal.genes$Gst2))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Gst2, get.gene.mode, target.region = gene.coords$Gst2[g, 2:3])) - 2
hmm.cnv.table[["Max_Gst2"]] <- apply(hmm.cnv.table[, names(focal.genes$Gst2)], 1, max)
cov.cnv.samples$Gst2 <- rownames(subset(hmm.cnv.table, Max_Gst2 > 0))

#Harb1
for (g in names(focal.genes$Harb1))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Harb1, get.gene.mode, target.region = gene.coords$Harb1[g, 2:3])) - 2
hmm.cnv.table[["Max_Harb1"]] <- apply(hmm.cnv.table[, names(focal.genes$Harb1)], 1, max)
cov.cnv.samples$Harb1 <- rownames(subset(hmm.cnv.table, Max_Harb1 > 0))

#Zinccarbo
for (g in names(focal.genes$Zinccarbo))
        hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$Zinccarbo, get.gene.mode, target.region = gene.coords$Zinccarbo[g, 2:3])) - 2
hmm.cnv.table[["Max_Zinccarbo"]] <- apply(hmm.cnv.table[, names(focal.genes$Zinccarbo)], 1, max)
cov.cnv.samples$Zinccarbo <- rownames(subset(hmm.cnv.table, Max_Zinccarbo > 0))


# Load up the 3RL
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


# Load up the data for X
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

# Set the diagnostic reads that will be used to detect each CNV
known.cnvs <- list(
  Ache = list(),
  Cyp6 = list(
    Dup1 = list(
      FA = matrix(c(8661000, 8714500, 8661400, 8714800), 2, 2),
      BP = data.frame(pos = c(8661135, 8714965), seq = c("AACGA", "ATGGC"))
    ),
    Dup2 = list(
      FA = matrix(c(8668500, 8674900, 8668800, 8675200), 2, 2),
      BP = data.frame(pos = c(8668667, 8675041), seq = c("GGTAG", "TTCAT"))
    ),
    Dup3 = list(
      FA = matrix(c(8669300, 8675600, 8669600, 8675900), 2, 2),
      BP = data.frame(pos = c(8669407, 8675890), seq = c("CACCA", "TGATG"))
    ),
    Dup4 = list(
      FA = matrix(c(8669700, 8675750, 8670000, 8676050), 2, 2),
      BP = data.frame(pos = c(8669896, 8675821), seq = c("CAAGT", "GGATG"))
    ),
    Dup5 = list(
      FA = matrix(c(8669700, 8686200, 8670000, 8686500), 2, 2),
      BP = data.frame(pos = c(8669807, 8686375), seq = c("TATTG", "GACAT"))
    ),
    Dup6 = list(
      FA = matrix(c(8669000, 8674400, 8669300, 8674700), 2, 2),
      BP = data.frame(pos = c(8669430, 8674523), seq = c("AATGT", "GCGCG"))
    ),
    Dup7 = list(
      FA = matrix(c(8671979, 8674100, 8672768, 8674658), 2, 2),
      BP = data.frame(pos = c(8672086, 8674715), seq = c("CCTCC", "GTTTT"))
    ),
    Dup8 = list(
      FA = matrix(c(8690250, 8696350, 8690550, 8696650), 2, 2),
      BP = data.frame(pos = c(8690385, 8696457), seq = c("CATGC", "TGAGC"))
    ),
    Dup9 = list(
      FA = matrix(c(8706331, 8709200, 8706700, 8709500), 2,2),
      BP = data.frame(pos = c(8706631, 8709386), seq = c("GAAAA", "CGTTT"))
    ),
    Insert_6.5 = list(
      FA = matrix(c(8403301, 8687323, 8404394, 8687823), 2,2),
      BP = data.frame(pos = c(8403989, 8410113), seq = c("TTTGC", "GTGCA"))
    )
  ),
  Esterase = list(),
  Gst1 = list(),
  Gst2 = list(
    Dup1 = list(
      FA = matrix(c(76398425, 76416209, 76401139, 76419386), 2, 2)#,
    #   BP = data.frame(pos = c(76402834, 76405258), seq = c("CTCAA", "TAGTA"))
    ),
    Dup2 = list(
      FA = matrix(c(76403487, 76409539, 76405933, 76411888), 2, 2)#,
    #   BP = data.frame(pos = c(76405273, 76411608), seq = c("AAGGT", "AGATA"))
    ),
    Dup3 = list(
      FA = matrix(c(76405506, 76408308, 76408765, 76409990), 2, 2),
      BP = data.frame(pos = c(76407222, 76409525), seq = c("GCTGC", "CGTTG"))
    ),
    Dup4 = list(
      FA = matrix(c(76662067, 76668449, 76663639, 76669851), 2, 2)#,
    #   BP = data.frame(pos = c(76667685, 76666005), seq = c("TGTTC", "TCAAA"))
    ),
    Dup5 = list(
      FA = matrix(c(76181149, 76195001, 76182080, 76196142), 2, 2)#,
    #   BP = data.frame(pos = c(76193011, 76194423), seq = c("TACAT", "AGTAT"))
    ),
    Dup6 = list(
      FA = matrix(c(76870049, 76875976, 76871359, 76876929), 2, 2)#,
    #   BP = data.frame(pos = c(76871233, 76872073), seq = c("GAGAA", "CATCG"))
    )
  ),
  Harb1 = list(),
  Zinccarbo = list(),
  GABA = list(),
  Carboxypep = list(),
  Cyp9 = list(),
  RDGA = list(),
  NADHCyp = list()
)

# In R < 4, characters are automatically converted to strings in data.frames, so fix that for all known.cnv
# entries here
for (region in names(known.cnvs)){
	for (Dup in names(known.cnvs[[region]])){
		if ('BP' %in% names(known.cnvs[[region]][[Dup]])){
			known.cnvs[[region]][[Dup]]$BP$seq <- as.character(known.cnvs[[region]][[Dup]]$BP$seq)
		}
	}
}

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

# Function to count the diagnostic reads for a given CNV. 
# di is a table of observed discordant and breakpoint reads for a given sample. 
count.reads.per.dup <- function(Dup.diagnostics, di){
	# There shouldn't be more than one of FA, SS, FM and XC
	if (sum(names(Dup.diagnostics) %in% c('FA', 'SS', 'FM', 'XC')) > 1)
		stop('There should not be more than one type of diagnostic read for a given Duplication.')
	num.reads <- setNames(numeric(7), c('FA', 'SS', 'FM', 'XCS', 'XCE', 'CSP', 'CEP'))
	if ('FA' %in% names(Dup.diagnostics)){
		num.reads['FA'] <- sum((di$FA$Position > Dup.diagnostics$FA[1,1]) & di$FA$Position < Dup.diagnostics$FA[1,2]
							 & (di$FA$Mate.position > Dup.diagnostics$FA[2,1]) & (di$FA$Mate.position < Dup.diagnostics$FA[2,2]))
	}
	if ('SS' %in% names(Dup.diagnostics)){
		num.reads['SS'] <- sum((di$SS$Position > Dup.diagnostics$SS[1,1]) & di$SS$Position < Dup.diagnostics$SS[1,2]
		                     & (di$SS$Mate.position > Dup.diagnostics$SS[2,1]) & (di$SS$Mate.position < Dup.diagnostics$SS[2,2]))
	}
	if ('FM' %in% names(Dup.diagnostics)){
		num.reads['FM'] <- sum((di$FM$Position > Dup.diagnostics$FM[1,1]) & di$FM$Position < Dup.diagnostics$FM[1,2]
							 & (di$FM$Mate.position > Dup.diagnostics$FM[2,1]) & (di$FM$Mate.position < Dup.diagnostics$FM[2,2]))
	}
	if ('XC' %in% names(Dup.diagnostics)){
		these.XC <- di$XC
		these.XC$Matechrom <- sub(':.*', '', these.XC$Mate.position)
		these.XC$Matepos <- as.integer(sub('.*:' , '', these.XC$Mate.position))
		num.reads['XCS'] <- sum((these.XC$Position > Dup.diagnostics$XC[1,1]) & (these.XC$Position < Dup.diagnostics$XC[1,2])
		                      & (these.XC$Matechrom == Dup.diagnostics$XC[2,3]) & (these.XC$Matepos > Dup.diagnostics$XC[2,1]) & (these.XC$Matepos < Dup.diagnostics$XC[2,2]), na.rm = T)
		num.reads['XCE'] <- sum((these.XC$Position > Dup.diagnostics$XC[3,1]) & (these.XC$Position < Dup.diagnostics$XC[3,2])
		                      & (these.XC$Matechrom == Dup.diagnostics$XC[4,3]) & (these.XC$Matepos > Dup.diagnostics$XC[4,1]) & (these.XC$Matepos < Dup.diagnostics$XC[4,2]), na.rm = T)
	}
	if ('BP' %in% names(Dup.diagnostics)){
		CEP.seq <- subset(di$BP$CEP, Position == Dup.diagnostics$BP$pos[1])$Clipped_sequence
		num.reads['CEP'] <- sum(substr(reverse(CEP.seq), 1, 5) == Dup.diagnostics$BP$seq[1])
		CSP.seq <- subset(di$BP$CSP, Position == Dup.diagnostics$BP$pos[2])$Clipped_sequence
		num.reads['CSP'] <- sum(substr(CSP.seq, 1, 5) == Dup.diagnostics$BP$seq[2])
	}
	num.reads
}

# Function to apply count.reads.per.dup to all CNVs for a given sample
count.diagnostic.reads <- function(diagnostic.reads.list, known.cnvs.list, all.breakpoint.positions, verbose = F){
	# Subset the BP reads so that they only include reads that match at least one of
	# the breakpoints. That will greatly reduce the list that needs searching for every
	# CNV. 
	diagnostic.reads.list$BP$CSP <- subset(diagnostic.reads.list$BP$CSP, Position %in% all.breakpoint.positions)
	diagnostic.reads.list$BP$CEP <- subset(diagnostic.reads.list$BP$CEP, Position %in% all.breakpoint.positions)
	num.reads <- lapply(known.cnvs.list, count.reads.per.dup, diagnostic.reads.list)
	if (verbose)
		print(num.reads)
	sapply(num.reads, sum)
}

# Function to apply counts.diagnostic.reads to all samples
count.diagnostic.reads.allsamples <- function(diagnostic.reads.allsamples, known.cnvs.list){
	if (length(known.cnvs.list) == 0)
		read.counts <- numeric()
	else{
		# Create a single object containing all of the breakpoint positions
		breakpoint.summary <- unlist(sapply(known.cnvs.list, function(x) x$BP$pos))
		# Re-wrote the following line function, because it failed when there was only one known CNV (outputs 
		# a table with the wrong orientation). 
		# t(sapply(diagnostic.reads.allsamples, count.diagnostic.reads, known.cnvs.list, breakpoint.summary))
		read.counts <- sapply(diagnostic.reads.allsamples, count.diagnostic.reads, known.cnvs.list, breakpoint.summary)
	}
	matrix(read.counts, nrow = length(diagnostic.reads.allsamples), ncol = length(known.cnvs.list), 
	       dimnames = list(names(diagnostic.reads.allsamples), names(known.cnvs.list)), byrow = T)
}

# Function to produce a table of all observed CNVs in all samples, given a list of coverage-based calls,
# a table of counts of diagnostic reads for each CNV in all samples, and a threshold number of diagnostic
# reads required to call a CNV
get.read.based.cnvs <- function(coverage.cnv.calls, diagnostic.read.counts.table, threshold.diagnostic.reads = 4){
	read.based <- cbind(diagnostic.read.counts.table >= threshold.diagnostic.reads)
	read.cnv.samples <- rownames(read.based)[apply(read.based, 1, any)]
	# Unlike in phase2, we use Dup0 to indicate that there is ANY coverage call, not just the samples that
	# have a coverage call but not a read-based call
	cbind(Dup0 = rownames(read.based) %in% coverage.cnv.calls, read.based)
}

# Set the folders in which to look for FA, SS, FM and XC reads. 
SSFA.folders <- setNames(paste(diagnostic.reads.folder, c('SSFA/2RL/Ache_region',
                                                          'SSFA/2RL/Cyp6_region',
                                                          'SSFA/2RL/Esterase_region',
                                                          'SSFA/2RL/Gst1_region',
                                                          'SSFA/2RL/Gst2_region',
                                                          'SSFA/2RL/Harb1_region',
                                                          'SSFA/2RL/Zinccarbo_region',
                                                          'SSFA/3RL/GABA_region',
                                                          'SSFA/3RL/Carboxypep_region',
                                                          'SSFA/X/Cyp9_region',
                                                          'SSFA/X/RDGA_region',
                                                          'SSFA/X/NADHCyp_region'),
                                                        sep = '/'),
                         c('Ache', 'Cyp6', 'Esterase', 'Gst1', 'Gst2', 'Harb1', 'Zinccarbo', 
						   'GABA', 'Carboxypep', 'Cyp9', 'RDGA', 'NADHCyp')
                        )

# Set the folders in which to look for breakpoint reads. 
breakpoints.folders <- gsub('SSFA', 'breakpoints', SSFA.folders)

cat('Loading discordant and breakpoint reads\n')
diagnostic.reads <- mcmapply(get.diagnostic.reads.allsamples, 
                             SSFA.folders,
                             breakpoints.folders, 
                             region.coords,
                             MoreArgs = list(sample.names = sample.names),
                             SIMPLIFY = F,
                             mc.preschedule = F, mc.cores = num.cores)

cat('Counting diagnostic reads\n')
diagnostic.read.counts <- mapply(count.diagnostic.reads.allsamples, 
                                 diagnostic.reads, known.cnvs)
cat('Calling read-based CNVs\n')
read.based.cnvs <- mapply(get.read.based.cnvs, cov.cnv.samples, diagnostic.read.counts, SIMPLIFY = F)

# Combine the separate CNV tables into a single table
full.cnv.table <- do.call(cbind, read.based.cnvs)
full.cnv.table <- cbind(rownames(full.cnv.table) %in% high.var.samples, full.cnv.table)
gene.cluster.names <- setNames(nm = names(diagnostic.read.counts))
# The column names for the full table will be with the gene cluster codes that we used in the phase 2 analysis.
colnames(full.cnv.table) <- c('High.var.sample', unlist(sapply(names(gene.cluster.names), function(x) paste(gene.cluster.names[x], colnames(read.based.cnvs[[x]]), sep = '_'))))
# Get the results in a different format. Here as a list where, for each CNV, we have a vector of samples that 
# carry it. 
cnv.sample.lists <- lapply(read.based.cnvs, function(m) apply(m, 2, function(x) rownames(m)[x]))

# Get a vector of samples that carry at least one CNV based on read calls.
# Dup0 is always the first column, so the -1 index removes Dup0
read.cnv.samples <- lapply(read.based.cnvs, function(x) rownames(x)[apply(x[, -1, drop = F], 1, any)])

# For each cluster, get a vector of samples that have been called as carrying a cnv by reads but not coverage,
# ignoring samples with high variance.
cov.based.negatives <- lapply(mapply(setdiff, read.cnv.samples, cov.cnv.samples), setdiff, high.var.samples)
# And vice-versa
read.based.negatives <- lapply(mapply(setdiff, cov.cnv.samples, read.cnv.samples), setdiff, high.var.samples)

dir.create('target_regions_analysis', showWarnings = FALSE)
write.table(full.cnv.table, file = 'target_regions_analysis/focal_region_CNV_table.csv', sep = '\t', col.names = NA, quote = F)
write.table(hmm.cnv.table, file = 'target_regions_analysis/HMM_gene_copy_number.csv', sep = '\t', col.names = NA, quote = F)

source(plotting.functions.file)
save.image('target_regions_analysis/target_regions_analysis.Rdata')

