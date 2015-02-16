#! /usr/bin/Rscript
# Plot the histogram of coverage by genomeCoverageBed
# plot-coverage.R [INPUT GENOMECOVERAGEBED] [WINDOW OVERALAP FLOAT] [REFERENCE]
# ./plot-coverage.R alignment.coverage 0.5 pXO1
library("ggplot2")

# Get arguements
args<-commandArgs(TRUE)
file <- args[1]
file_without_ext <- sub("^([^.]*).*", "\\1", file)
window_overlap <- as.double(args[2])
window_sizes <- c(100, 500, 1000, 2000, 5000, 10000)
reference <- args[3]

# Read input coverage file
coverage <- read.table(file, sep="\t", header=FALSE)
colnames(coverage) <- c("reference", "position", "coverage")
n <- length(coverage$coverage)

# Produce plots
output <- paste(file_without_ext, ".pdf", sep="")
pdf(output, onefile=TRUE, width=14, height=4)
for (j in 1:length(window_sizes)) {
    window_size <- window_sizes[j]
    windows <- seq(1, n-window_size, by=(window_size * window_overlap))
    print(paste("Creating plot for window size:", window_size))
    mean_per_window = sapply(1:length(windows), function(i) {
        mean(coverage$coverage[windows[i]:(windows[i]+(window_size-1))])
    })
    title <- paste("Coverage accross", reference,"for input",
                   basename(file_without_ext), "using sliding window of size:",
                   window_size, "bp (overlap", (window_size * window_overlap),
                   "bp)")
    p <- ggplot(data.frame(x=windows, y=mean_per_window), aes(x=x, y=y)) +
        geom_bar(stat="identity") +
        ggtitle(title) +
        xlab("Position (bp)") +
        ylab("Coverage") +
        theme_bw()

    # Save plot
    print(p)
}
dev.off()
