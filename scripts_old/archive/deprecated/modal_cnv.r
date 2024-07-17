
arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

chroms = c('2L', '2R', '3L', '3R', 'X')

threshold.variance <- 0.2

# Get the sample manifest
#sample.list <- arg.values[1]
sample.list <- '~/personal/phase3_data/phase3_manifest_arabiensis.txt'
sample.names <- setNames(nm = read.table(sample.list, stringsAsFactors = F)[[1]])

# Get the table of gene coordinates
#gene.regions.file <- arg.values[2]
gene.regions.file <- '~/personal/phase3_data/tables/gene_regions.csv'
all.gene.coordinates <- read.table(gene.regions.file, sep = '\t', header = T, row.names = 1)
gene.coordinates.list <- lapply(split(all.gene.coordinates, all.gene.coordinates$Chrom)[chroms], function(x) x[, c('start', 'end')])

# Get the metadata
#meta.file <- arg.values[3]
meta.file <- '~/personal/phase3_data/phase3.samples.meta.csv'
meta <- read.table(meta.file, sep = '\t', header = T, row.names = 1, quote = '', comment.char = '')[sample.names, ]
expected.copy.number.on.X <- c(2, 1)[(meta$sex_call == 'M') + 1]

# Get the species name (only used to name the output file)
species <- arg.values[4]

# Get the coverage variance data
#cov.var.file <- arg.values[5]
cov.var.file <- '~/personal/phase3_cnv/coverage/coverage_variance_masked_09_05_all.csv'
cov.var <- read.table(cov.var.file, header = T, sep = '\t', row.names = 1)
high.var.samples <- intersect(sample.names, rownames(cov.var)[cov.var$autosomes > threshold.variance])

# Get the folder in which the HMM output are stored
#coverage.folder <- arg.values[6]

# Get the output folder
output.folder <- arg.values[7]
dir.create(output.folder, showWarnings = FALSE)

# Function to load up each file and store the results in a list
load.hmm.files <- function(sample.names, chrom){
	folder <- paste(coverage.folder, chrom, 'HMM_output', sep = '/')
	all.files <- setNames(paste(folder, '/', sample.names, '_', chrom, '_HMM.csv', sep = ''), sample.names)
	output.table <- as.data.frame(lapply(all.files, function(fn) read.table(fn, header = T, row.names = 1)[, 'CNV', drop = F]))
	colnames(output.table) <- sample.names
	output.table$Position <- as.numeric(rownames(output.table))
	output.table
}

# Functions that will calculate the mode of the HMM coverage in a given region
calculate.mode <- function(copy.number.vector){
	# We calculate the mode. Where there is more than one mode, we take the smallest one
	un <- unique(copy.number.vector)
	tab <- tabulate(match(copy.number.vector, un))
	copy.number.mode <- min(un[tab == max(tab)])
	copy.number.mode
}

get.gene.mode <- function(target.region, hmm.data, window.size = 300){
	target.region <- as.numeric(target.region)
	hmm.in.region <- subset(hmm.data, Position >= (target.region[1] - window.size) & Position < target.region[2])[, sample.names]
	gene.mode <- apply(hmm.in.region, 2, calculate.mode)
	gene.mode
}

modal.copy.number.change <- function(chrom, expected.copy.number){
	cat(chrom, '\n')
	cat('\tLoading hmm output files.\n')
	full.table <- load.hmm.files(sample.names, chrom)
	cat('\tCalculating modal copy number change.\n')
	modal.copy.number.change <- apply(gene.coordinates.list[[chrom]], 1, get.gene.mode, hmm.data = full.table) - expected.copy.number
}

modal.copy.number <- do.call(cbind, mapply(modal.copy.number.change, chroms, list(2, 2, 2, 2, expected.copy.number.on.X)))

modal.copy.number <- modal.copy.number[, order(colnames(modal.copy.number))]

# Write that table to file
output_fn = paste(output_folder, '/modal_copy_number_', species, '.csv', sep = '')
write.table(modal.copy.number, file = output_fn, sep = '\t', col.names = NA, quote = F)

