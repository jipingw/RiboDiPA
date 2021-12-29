\name{exonTrack}
\alias{exonTrack}
\title{
Visualization: generating tracks for \code{igvR} (per exon)
}
\description{
This function outputs a list of \code{GRanges} format objects of 
ribosome P-site footprints per exon, as well as corresponding test results. 
It can be used for \code{igvR} visualization.
}
\usage{
exonTrack(data, gene)
}
\arguments{
\item{data}{Data object from \code{diffPatterbTestExon} function or wrapper
function \code{RiboDiPA}.
}
\item{gene}{One gene for visualization at a time, since different genes have 
different number of transcripts.
}
}
\details{
R package \code{igvR} is not implemented in \code{RiboDiPA}. Users need install
\code{igvR} through Bioconductor or relevant source package. A simple
illustration example of how to use it for \code{igvR} visualization is given
below.
}
\value{
A list of lists. Each element is a list of \code{GRanges} objects representing
replicates. Each second level list element is P-site footprint count per exon 
with differential pattern test results in a transcript.
}
\examples{
data(result.exon)
tracks.exon <- exonTrack(data = result.exon, gene = "tY(GUA)D")

# library(igvR)
# igv <- igvR()
# setBrowserWindowTitle(igv, "ribosome footprint per exon example")
# setGenome(igv, "saccer3")

# for(track.name in names(tracks.exon)){
#     track.rep <- tracks.exon[[track.name]]
#     for(tx.name in names(track.rep)){
#         track.tx <- tracks.exon[[track.name]][[tx.name]]
#         track <- GRangesQuantitativeTrack(trackName = 
#             paste(track.name, tx.name), track.tx[,1], color = track.name)
#         displayTrack(igv, track)
#     }
# }
}
