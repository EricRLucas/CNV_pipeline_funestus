gene.colours <- c('black', 'purple', 'orange', 'magenta', 'brown', 'green', 'violet', 'red', 'blue')

# Create a function that can be used to plot discordant read pairs
add.diagnostics <- function(coordinates, this.col = 'red', yrange = c(0,1)){
	coordinates <- as.matrix(coordinates)
	n <- nrow(coordinates)
	co <- ncol(coordinates)
	jitter.values <- seq(yrange[1], yrange[2], length.out = n)
	points(coordinates, matrix(jitter.values, nrow = n, ncol = co), col = this.col)
	if (co == 2)
		segments(coordinates[,1], jitter.values, coordinates[,2], jitter.values, col = this.col)
}


plot.Cyp6 <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Cyp6[1], end.pos = plotting.ranges$Cyp6[2]){
	start.index <- which(compact.hmm$Cyp6[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Cyp6[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Cyp6[[this.sample]]$Position[start.index : end.index], compact.hmm$Cyp6[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Cyp6[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Cyp6[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Cyp6[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Cyp6[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Cyp6[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Cyp6[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Cyp6)]
	abline(v = unlist(gene.coords$Cyp6[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Cyp6[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Cyp6), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Cyp6 <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Cyp6[1], end.pos = plotting.ranges$Cyp6[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Cyp6(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Ache <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Ache[1], end.pos = plotting.ranges$Ache[2]){
	start.index <- which(compact.hmm$Ache[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Ache[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Ache[[this.sample]]$Position[start.index : end.index], compact.hmm$Ache[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Ache[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Ache[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Ache[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Ache[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Ache[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Ache[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Ache)]
	abline(v = unlist(gene.coords$Ache[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Ache[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Ache), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Ache <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Ache[1], end.pos = plotting.ranges$Ache[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Ache(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Esterase <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Esterase[1], end.pos = plotting.ranges$Esterase[2]){
	start.index <- which(compact.hmm$Esterase[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Esterase[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Esterase[[this.sample]]$Position[start.index : end.index], compact.hmm$Esterase[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Esterase[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Esterase[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Esterase[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Esterase[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Esterase[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Esterase[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Esterase)]
	abline(v = unlist(gene.coords$Esterase[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Esterase[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Esterase), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Esterase <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Esterase[1], end.pos = plotting.ranges$Esterase[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Esterase(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Gst1 <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Gst1[1], end.pos = plotting.ranges$Gst1[2]){
	start.index <- which(compact.hmm$Gst1[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Gst1[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Gst1[[this.sample]]$Position[start.index : end.index], compact.hmm$Gst1[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Gst1[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Gst1[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Gst1[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Gst1[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Gst1[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Gst1[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Gst1)]
	abline(v = unlist(gene.coords$Gst1[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Gst1[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Gst1), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Gst1 <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Gst1[1], end.pos = plotting.ranges$Gst1[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Gst1(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Harb1 <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Harb1[1], end.pos = plotting.ranges$Harb1[2]){
	start.index <- which(compact.hmm$Harb1[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Harb1[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Harb1[[this.sample]]$Position[start.index : end.index], compact.hmm$Harb1[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Harb1[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Harb1[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Harb1[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Harb1[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Harb1[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Harb1[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Harb1)]
	abline(v = unlist(gene.coords$Harb1[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Harb1[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Harb1), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Harb1 <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Harb1[1], end.pos = plotting.ranges$Harb1[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Harb1(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Gst2 <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Gst2[1], end.pos = plotting.ranges$Gst2[2]){
	start.index <- which(compact.hmm$Gst2[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Gst2[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Gst2[[this.sample]]$Position[start.index : end.index], compact.hmm$Gst2[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Gst2[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Gst2[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Gst2[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Gst2[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Gst2[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Gst2[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Gst2)]
	abline(v = unlist(gene.coords$Gst2[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Gst2[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Gst2), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Gst2 <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Gst2[1], end.pos = plotting.ranges$Gst2[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Gst2(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Zinccarbo <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Zinccarbo[1], end.pos = plotting.ranges$Zinccarbo[2]){
	start.index <- which(compact.hmm$Zinccarbo[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Zinccarbo[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Zinccarbo[[this.sample]]$Position[start.index : end.index], compact.hmm$Zinccarbo[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Zinccarbo[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Zinccarbo[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Zinccarbo[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Zinccarbo[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Zinccarbo[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Zinccarbo[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Zinccarbo)]
	abline(v = unlist(gene.coords$Zinccarbo[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Zinccarbo[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Zinccarbo), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Zinccarbo <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Zinccarbo[1], end.pos = plotting.ranges$Zinccarbo[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Zinccarbo(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.GABA <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$GABA[1], end.pos = plotting.ranges$GABA[2]){
	start.index <- which(compact.hmm$GABA[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$GABA[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$GABA[[this.sample]]$Position[start.index : end.index], compact.hmm$GABA[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$GABA[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$GABA[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$GABA[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$GABA[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$GABA[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$GABA[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$GABA)]
	abline(v = unlist(gene.coords$GABA[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$GABA[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$GABA), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.GABA <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$GABA[1], end.pos = plotting.ranges$GABA[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.GABA(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Carboxypep <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Carboxypep[1], end.pos = plotting.ranges$Carboxypep[2]){
	start.index <- which(compact.hmm$Carboxypep[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Carboxypep[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Carboxypep[[this.sample]]$Position[start.index : end.index], compact.hmm$Carboxypep[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Carboxypep[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Carboxypep[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Carboxypep[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Carboxypep[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Carboxypep[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Carboxypep[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Carboxypep)]
	abline(v = unlist(gene.coords$Carboxypep[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Carboxypep[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Carboxypep), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Carboxypep <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Carboxypep[1], end.pos = plotting.ranges$Carboxypep[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Carboxypep(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.Cyp9 <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Cyp9[1], end.pos = plotting.ranges$Cyp9[2]){
	start.index <- which(compact.hmm$Cyp9[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$Cyp9[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$Cyp9[[this.sample]]$Position[start.index : end.index], compact.hmm$Cyp9[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$Cyp9[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$Cyp9[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$Cyp9[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$Cyp9[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$Cyp9[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$Cyp9[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$Cyp9)]
	abline(v = unlist(gene.coords$Cyp9[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$Cyp9[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$Cyp9), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.Cyp9 <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$Cyp9[1], end.pos = plotting.ranges$Cyp9[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.Cyp9(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.RDGA <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$RDGA[1], end.pos = plotting.ranges$RDGA[2]){
	start.index <- which(compact.hmm$RDGA[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$RDGA[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$RDGA[[this.sample]]$Position[start.index : end.index], compact.hmm$RDGA[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$RDGA[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$RDGA[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$RDGA[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$RDGA[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$RDGA[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$RDGA[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$RDGA)]
	abline(v = unlist(gene.coords$RDGA[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$RDGA[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$RDGA), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.RDGA <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$RDGA[1], end.pos = plotting.ranges$RDGA[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.RDGA(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


plot.NADHCyp <- function(this.sample, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$NADHCyp[1], end.pos = plotting.ranges$NADHCyp[2]){
	start.index <- which(compact.hmm$NADHCyp[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$NADHCyp[[this.sample]]$Position <= end.pos), 1)
	if (this.sample %in% high.var.samples)
		par(bg = 'beige')
	else
		par(bg = 'white')
	plot(compact.hmm$NADHCyp[[this.sample]]$Position[start.index : end.index], compact.hmm$NADHCyp[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$NADHCyp[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$NADHCyp[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$NADHCyp[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$NADHCyp[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$NADHCyp[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$NADHCyp[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$NADHCyp)]
	abline(v = unlist(gene.coords$NADHCyp[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$NADHCyp[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$NADHCyp), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.NADHCyp <- function(list.of.samples, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$NADHCyp[1], end.pos = plotting.ranges$NADHCyp[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot.NADHCyp(this.sample, diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


simple.plot <- function(sample.names, gene.cluster, smoothing = 5, maxy = 12, ...){
	these.coords <- plotting.ranges[[gene.cluster]]
	plot(c(min(these.coords), max(these.coords)), c(0,maxy), type = 'n', bty = 'n', xaxt = 'n', xlab = '', ylab = 'Normalised coverage')
	axis(1, lwd = 0, mgp = c(0,0.25,2))
	mtext(paste('Position on chromosome', unique(gene.coords[[gene.cluster]]$Chrom)), 1, 1.75)
	these.gene.colours <- gene.colours[1:nrow(gene.coords[[gene.cluster]])]
	abline(v = unlist(gene.coords[[gene.cluster]][, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords[[gene.cluster]][, c('start', 'end')], 1, mean), 9, rownames(gene.coords[[gene.cluster]]), srt=90, adj = 0, col = these.gene.colours)
	smoothed.line <- function(s){
		these.counts <- subset(compact.hmm[[gene.cluster]][[s]], Position >= these.coords[1] & Position <= these.coords[2])
		# If needed, create a table of smoothed data here
		if (smoothing > 1){
			if (smoothing > nrow(these.counts))
				stop('Fail. Smoothing window size is larger than the number of points to smooth over.')
			these.counts <- t(sapply(1:(nrow(these.counts) - smoothing + 1), function(j) apply(these.counts[j:(j + smoothing - 1), c('Position', 'Normalised_coverage')], 2, mean)))
		}
		lines(these.counts[, 'Position'], these.counts[, 'Normalised_coverage'])
	}
	empty <- sapply(sample.names, smoothed.line)
}
