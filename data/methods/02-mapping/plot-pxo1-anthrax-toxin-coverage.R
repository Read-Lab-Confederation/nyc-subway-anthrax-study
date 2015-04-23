#! /usr/bin/Rscript
# Plot the histogram of mean coverage (1000bp windows, 50% overlap) for the
# complete pXO1 plasmid. Then creates subplots of the actual read coverage of
# PXO1 lethal genes (cya, lef, pagA, pagR).
#
# plot-pxo1-anthrax-toxin-coverage.R [INPUT GENOMECOVERAGEBED] [GFF3 ANNOTATION]
# ./plot-pxo1-anthrax-toxin-coverage.R pXO1.coverage pXO1.gff3
library("Gviz")
library("rtracklayer")
options(ucscChromosomeNames = FALSE)

get_gene_tracks <- function(gene, name, coverage, track_length, max_coverage) {
    flank <- track_length - width(gene)
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
pxo1_gff <- import.gff3(
    args[2],
    genome= "PXO1",
    asRangedData = FALSE,
    feature.type = c("CDS")
)

# Whole pXO1 Coverage Track (1000bp sliding window 500 bp overlap)
window_size <- 1000
window_overlap <- 0.5
windows <- seq(1, n-window_size, by=(window_size * window_overlap))
mean_per_window = sapply(1:length(windows), function(i) {
    mean(coverage$coverage[windows[i]:(windows[i]+(window_size-1))])
})

# Get genes
cya <- pxo1_gff[elementMetadata(pxo1_gff)[,"Name"]=="cya"]
lef <- pxo1_gff[elementMetadata(pxo1_gff)[,"Name"]=="lef"]
paga <- pxo1_gff[elementMetadata(pxo1_gff)[,"Name"]=="pagA"]
pagr <- pxo1_gff[elementMetadata(pxo1_gff)[,"Name"]=="pagR"]

max_coverage <- max(
    mean_per_window,
    coverage$coverage[start(cya):end(cya)],
    coverage$coverage[start(lef):end(lef)],
    coverage$coverage[start(paga):end(paga)],
    coverage$coverage[start(pagr):end(pagr)]
)

# Get genes tracks (axis, coverage, name)
scale_width <- max(width(lef), width(cya),width(paga), width(pagr))

cya <- get_gene_tracks(cya, "cya", coverage, scale_width + 200, max_coverage)
lef <- get_gene_tracks(lef, "lef", coverage, scale_width + 200, max_coverage)
paga <- get_gene_tracks(paga, "pagA", coverage, scale_width + 200, max_coverage)
pagr <- get_gene_tracks(pagr, "pagR", coverage, scale_width + 200, max_coverage)

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
        start = c(cya$start, lef$start, paga$start, pagr$start),
        end = c(cya$end, lef$end, paga$end, pagr$end),
        names = c("cya", "lef", "pagA", "pagR")
))

# Produce plots
# PDF
output <- paste(file_without_ext, "-anthrax-toxin.pdf", sep="")
pdf(output, onefile=TRUE, width=12, height=16)
grid.newpage()
pushViewport(viewport(height=0.20, y=1, just="top"))
plotTracks(list(complete_axis, complete_coverage), add=TRUE, cex.axis = 1,
           fill.range = c("red", "blue", "darkgreen", "darkorange"))
popViewport()
pushViewport(viewport(height=0.80, y=0, just="bottom", layout = grid.layout(4,1)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(list(lef$axis, lef$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("red"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plotTracks(list(pagr$axis, pagr$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("blue"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
plotTracks(list(paga$axis, paga$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("darkgreen"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 4))
plotTracks(list(cya$axis, cya$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("darkorange"))
popViewport(1)
popViewport()
dev.off()

# PNG
output_svg <- paste(file_without_ext, "-anthrax-toxin.svg", sep="")
svg(output_svg, width=12, height=16)
grid.newpage()
pushViewport(viewport(height=0.20, y=1, just="top"))
plotTracks(list(complete_axis, complete_coverage), add=TRUE, cex.axis = 1,
           fill.range = c("red", "blue", "darkgreen", "darkorange"))
popViewport()
pushViewport(viewport(height=0.80, y=0, just="bottom", layout = grid.layout(4,1)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(list(lef$axis, lef$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("red"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plotTracks(list(pagr$axis, pagr$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("blue"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
plotTracks(list(paga$axis, paga$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("darkgreen"))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 4))
plotTracks(list(cya$axis, cya$coverage), showId = TRUE,  cex.axis = 1,
           add=TRUE, fill.range = c("darkorange"))
popViewport(1)
popViewport()
dev.off()
