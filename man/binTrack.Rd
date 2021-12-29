\name{binTrack}
\alias{binTrack}
\title{
Visualization: generating tracks for \code{igvR} (per bin)
}
\description{
This function outputs a list of \code{GRanges} format objects of binned
ribosome P-site footprints, as well as test results. It can be used for
\code{igvR} visualization.
}
\usage{
binTrack(data, exon.anno)
}
\arguments{
\item{data}{Data object from \code{diffPatterbTest} function or wrapper
function \code{RiboDiPA}.
}
\item{exon.anno}{
\code{exons} value from \code{psiteMapping} function. Genomic start and 
end coordinates, and other gene information is included. See 
\code{psiteMapping} for details.
}
}
\details{
R package \code{igvR} is not implemented in \code{RiboDiPA}. Users need install
\code{igvR} through Bioconductor or relevant source package. A simple
illustration example of how to use it for \code{igvR} visualization is given
below.
}
\value{
A list of \code{GRanges} format objects. Each element is a \code{GRanges} 
object of the binned P-site footprint tracks with differential pattern test
results.
}
\examples{
data(result.pst)
data(data.psite)
tracks.bin <- binTrack(data = result.pst, 
    exon.anno = data.psite$exons)

# library(igvR)
# thred <- 0.05
# igv <- igvR()
# setBrowserWindowTitle(igv, "ribosome binned track example")
# setGenome(igv, "saccer3")

# for(track.name in names(tracks.bin)){
#     track.rep <- tracks.bin[[track.name]]
#     resize(track.rep, width(track.rep) + 1)
#     track <- GRangesQuantitativeTrack(trackName = paste(track.name,"binned"),
#         track.rep[,1], color = "black")
#     displayTrack(igv, track)
# }
# track.rep2 <- tracks.bin[[1]]
# sig.bin <- (values(track.rep2)[,5] <= thred)
# log10.padj <- - log10(values(track.rep2)[,5]) 
# mcols(track.rep2) <- data.frame(log10padj = log10.padj)
# track.rep2 <- track.rep2[which(sig.bin),]
# track <- GRangesQuantitativeTrack(trackName = "- log 10 of padj", 
#     track.rep2, color = "red", trackHeight = 40)
# displayTrack(igv, track)
}