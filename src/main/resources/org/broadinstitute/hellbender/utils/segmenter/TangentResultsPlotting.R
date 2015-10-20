# Useful for debugging
#options(error = quote({dump.frames(dumpto = "CBS_dump", to.file = TRUE); q()}))

option_list <- list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--output_file", "-output_file"), dest="output_file", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"))

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
output_file=opt[["output_file"]]
log_input=as.logical(opt[["log2_input"]])

# Use a function for debugging purposes
create_tangent_plots_file = function(sample_name, tn_file, output_file, log_input) {
	# Read in file and extract needed data
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	contig = tn[,"contig"]
	pos = tn[,"stop"]
	# Ability to switch between copy-ratio and log2 copy-ratio
	if (log_input) {
	    dat = tn[,sample_name]
	} else {
	    dat = log2(tn[,sample_name])
	}

}

FullGenPlots <- function(seg.dat, plots.dir, SID, bound_name=NULL, y.max=5) {
  cat("\nPlotting Full Genome...\n")
  plot.fn = file.path( plots.dir, paste(SID, "_FullGenome", bound_name, ".png", sep=""))
  png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
  chr.dat = SetUpPlot("Total copy ratio", 0, y.max, "", F)
  TotalCopyRatio(seg.dat, chr.dat)
  dev.off()
}

## Chromosome Background for data
SetUpPlot <- function(y.lab, y.min, y.max, x.lab, lab.chr){
  ## Via http://genome.ucsc.edu/cgi-bin/hgTables
  ## HG19 Chr Lens
  chr.lens <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)

  ## HG19 Cent Pos left and right bounds
  cent.posS <- c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679,
	39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782,
	26369569, 11288129, 13000000)
  cent.posE <- c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679,
	42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782,
	29369569, 14288129, 16000000)

  ###Add sex chromosomes###
  ## Add X
  #chr.lens = c(chr.lens, 155270560)
  #cent.posS = c(cent.posS, 58632012)
  #cent.posE = c(cent.posE, 61632012)

  ## Add Y - DISABLED
  #chr.lens = c(chr.lens, 59373566)
  #cent.posS = c(cent.posS, 10104553)
  #cent.posE = c(cent.posE, 13104553)

  chr.w <- chr.lens / sum(chr.lens)
  suppressWarnings( par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) )
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
  ww <- as.vector(rbind(chr.w, chr.w)) / 2
  chr.mids <- cumsum(ww)[(1:length(ww) - 1) %% 2 == 0]
  if(lab.chr){
  	lab.vals <- (c(1:length(chr.w)))
  	odd.ix <- lab.vals %% 2 == 1
  	mtext(text = lab.vals[odd.ix], side = 1, line = -0.45, at = chr.mids[odd.ix],
    		las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  	mtext(text = lab.vals[!odd.ix], side = 1, line = 0, at = chr.mids[!odd.ix],
      	las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  }
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  cent.posS <- cent.posS/sum(chr.lens) + chr.offsets[1:(length(chr.offsets) - 1)]
  cent.posE <- cent.posE/sum(chr.lens) + chr.offsets[1:(length(chr.offsets) - 1)]
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col <- ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1],
    		ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(cent.posS[i], 2), lty = 3, lwd = 0.5)
    lines(y = c(y.min, y.max), x = rep(cent.posE[i], 2), lty = 3, lwd = 0.5)
  }
  chr.dat = list(chr.w, chr.lens, chr.offsets)
  return( chr.dat )
}

## Plotting Copy Ratio
TotalCopyRatio <- function(seg.dat, chr.dat){
  chr.w = chr.dat[[1]]
  chr.lens = chr.dat[[2]]
  chr.offsets = chr.dat[[3]]
  for( s in 1:length(seg.dat[["as.res"]][["h.seg.dat"]][["h.capseg.d"]])) {
    Theta = seg.dat[["capture.em.fit"]][["Theta"]]
    delta.tau = seg.dat[["capture.em.fit"]][["delta.tau"]]
    atten.tau1 = AffyAtten(delta.tau[s, 2], Theta[["at.capseg"]])
    chr <- seg.dat$as.res$h.seg.dat$h.capseg.annot[[s]][["chr"]]
    short <- seg.dat[["as.res"]][['h.seg.dat']]
    seg.crds <- short[["h.capseg.annot"]][[s]][["pos"]]
    genome.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]

    colors=c("coral", "dodgerblue")

    points(genome.crds,short[["h.capseg.d"]][[s]][1,], col=colors[s%%2+1], pch=16, cex=0.2)
    lines(range(genome.crds), rep(atten.tau1, 2), col="black", lwd=1, lty=2)
  }
}

create_tangent_plots_file(sample_name, tn_file, output_file, log_input)
