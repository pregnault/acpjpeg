\name{acpjpeg-package}
\alias{acpjpeg-package}
\alias{acpjpeg}
\docType{package}
\title{
\packageTitle{acpjpeg}
}
\description{
\packageDescription{acpjpeg}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{acpjpeg}
\packageIndices{acpjpeg}
}
\author{
\packageAuthor{acpjpeg}

Maintainer: \packageMaintainer{acpjpeg}
}
\examples{
data(lenna)
jpeg::writeJPEG(lenna, target='lenna.jpeg')
pcajpeg(file="lenna.jpeg")
}
