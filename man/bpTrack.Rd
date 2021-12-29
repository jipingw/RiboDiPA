\name{bpTrack}
\alias{bpTrack}
\title{
Visualization: generating tracks for \code{igvR} (per bp)
}
\description{
This function outputs a list of \code{GRanges} format objects of
ribosome P-site footprints per bp. It can be used for \code{igvR}
visualization.
}
\usage{
bpTrack(data, names.rep = NULL, genes.list)
}
\arguments{
\item{data}{
Data object from \code{psiteMapping} function or wrapper function
\code{RiboDiPA}.
}
\item{names.rep}{
Customized names of the replicates. Default value uses the column names
of \code{data$counts}.
}
\item{genes.list}{
A list of genes for visualization. Default value uses all genes listed in
\code{data$coverage}
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
object of the P-site footprint tracks of selected genes.
}
\examples{
data(data.psite)
names.rep <- c("mutant1", "mutant2", "wildtype1", "wildtype2")
tracks.bp <- bpTrack(data = data.psite, names.rep = names.rep, 
    genes.list = c("YDR050C", "YDR062W", "YDR064W"))

# library(igvR)
# igv <- igvR()
# setBrowserWindowTitle(igv, "ribosome footprint track example")
# setGenome(igv, "saccer3")
# for(track.name in names.rep){
#     track.rep <- tracks.bp[[track.name]]
#     track <- GRangesQuantitativeTrack(trackName = paste(track.name, "bp"),
#         track.rep[,1], color = "green")
#     displayTrack(igv, track)
# }
}
