#! /usr/bin/Rscript
# Plot the histogram of mean coverage (1000bp windows, 50% overlap) for the
# complete pMT1 plasmid. Then create subplots of the actual read coverage of
# murine toxin genes.
#
# plot-pmt1-murine-toxin-coverage.R [INPUT GENOMECOVERAGEBED] [GFF3 ANNOTATION]
# ./plot-pmt1-murine-toxin-coverage.R pmt1.coverage pmt1.gff3
library("Gviz")
library("rtracklayer")
options(ucscChromosomeNames = FALSE)

get_gene_tracks <- function(gene, name, coverage, flank, max_coverage) {
    window <- seq(start(gene) - flank, end(gene) + flank)
    coverage = sapply(1:length(window), function(i) {
        as.double(coverage$coverage[window[i]])
    })

    coverage_track <- DataTrack(
        start = window,
        width = 1,
        name = " ",
        data = coverage,
        genome = "*",
        chromosome = "*",
        type="hist",
        ylim = c(0, max_coverage)
    )

    axis_track <- GenomeAxisTrack(
        labelPos = "above",
        range = IRanges(
            start = c(start(gene)), end = c(end(gene)), names = c(name)
    ))

    return(list(coverage=coverage_track, axis=axis_track, name=name,
                start=start(gene), end=end(gene)))
}

# Get arguements
args<-commandArgs(TRUE)
file = args[1]
file_without_ext <- gsub(".coverage.gz$", "", file)
coverage <- read.table(gzfile(file), sep="\t", header=FALSE)
colnames(coverage) <- c("reference", "position", "coverage")
n <- length(coverage$coverage)


# Read GFF3
pmt1_gff <- import.gff3(
    args[2],
    genome= "pMT1",
    asRangedData = FALSE,
    feature.type = c("CDS")
)

# Whole pMT1 Coverage Track (1000bp sliding window 500 bp overlap)
window_size <- 1000
window_overlap <- 0.5
windows <- seq(1, n-window_size, by=(window_size * window_overlap))
mean_per_window = sapply(1:length(windows), function(i) {
    mean(coverage$coverage[windows[i]:(windows[i]+(window_size-1))])
})

# Get ymt (in GenBank it is labled as gene=YPMT1.74, product=murine toxin)
ymt <- pmt1_gff[elementMetadata(pmt1_gff)[,"Name"]=="YPMT1.74"]

max_coverage <- max(
    mean_per_window,
    coverage$coverage[start(ymt):end(ymt)]
)

# Get genes tracks (axis, coverage, name)
ymt <- get_gene_tracks(ymt, "ymt", coverage, 200, max_coverage)

# Complete Plasmid Axis and Coverage
complete_coverage <- DataTrack(
    start = windows,
    width = 1,
    name = "Coverage",
    data = mean_per_window,
    genome = "*",
    chromosome = "*",
    type = "hist",
    ylim = c(0, max_coverage)
)

complete_axis <- GenomeAxisTrack(
    labelPos = "above",
    range = IRanges(
        start = c(ymt$start),
        end = c(ymt$end),
        names = c("ymt")
))

# Produce plots
# PDF
output <- paste(file_without_ext, "-murine-toxin.pdf", sep="")
pdf(output, onefile=TRUE, width=12, height=8)
grid.newpage()
pushViewport(viewport(height=0.50, y=1, just="top"))
plotTracks(list(complete_axis, complete_coverage), add=TRUE,
           fill.range = c("red"))
popViewport()
pushViewport(viewport(height=0.50, y=0, just="bottom"))
plotTracks(list(ymt$axis, ymt$coverage), showId = TRUE,
           add=TRUE, fill.range = c("red"))
popViewport()
dev.off()

# SVG
output_svg <- paste(file_without_ext, "-murine-toxin.svg", sep="")
svg(output_svg, width=12, height=8)
grid.newpage()
pushViewport(viewport(height=0.50, y=1, just="top"))
plotTracks(list(complete_axis, complete_coverage), add=TRUE,
           fill.range = c("red"))
popViewport()
pushViewport(viewport(height=0.50, y=0, just="bottom"))
plotTracks(list(ymt$axis, ymt$coverage), showId = TRUE,
           add=TRUE, fill.range = c("red"))
popViewport()
dev.off()

