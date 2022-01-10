#### Heading ####

## File: ACPJPEG.R
## Desc: Performs PCA on RGB channels of a JPEG image and returns the 
##       projections of the image on the first and two first components.
## Date: 13/10/15
## R version: 3.2.2
## Required packages: jpeg (>= 0.1-8), FactoMineR (>= 1.31.3), Rcpp (>= 0.12.1)
## Suggested packages: foreach (>= 1.4.2), microbenchmark (>= 1.4-2)
## philippe.regnault(at)univ-reims.fr


#### Coputing mean and standard deviation of RGBchannels ####

## Version degueux du calcul de la moyenne et de l'ecart-type des canaux de couleurs
CalcMS1 <- function(imgAsMat) {
  R <- imgAsMat[, 1]
  G <- imgAsMat[, 2]
  B <- imgAsMat[, 3]
  MR <- 0
  MG <- 0
  MB <- 0
  for (i in 1:length(R)) {
    MR <- MR + R[i]
    MG <- MG + G[i]
    MB <- MB + B[i]
  }
  M <- c(MR/length(R), MG/length(R), MB/length(R))
  
  SR <- 0
  SG <- 0
  SB <- 0
  for (i in 1:length(R)) {
    SR <- SR + (R[i] - M[1]) * (R[i] - M[1])
    SG <- SG + (G[i] - M[2]) * (G[i] - M[2])
    SB <- SB + (B[i] - M[3]) * (B[i] - M[3])
  }
  S <- sqrt(c(SR/length(R), SG/length(G), SB/length(B)))
  return(cbind(M,S))
}

## Version vectorisee en utilisant apply
CalcMS2 <- function(imgAsMat) {
  M <- apply(imgAsMat, 2, mean)
  S <- apply(imgAsMat, 2, sd)
  return(cbind(M,S))
}

## Version parallelisee avec mclapply (pour linux uniquement)
CalcMS3 <- function(imgAsLst) {
  M <- unlist(mclapply(imgAsLst, mean, mc.preschedule=FALSE, mc.cores=4))
  S <- unlist(mclapply(imgAsLst, sd, mc.preschedule=FALSE, mc.cores=4))
  return(cbind(M,S))
}

## Version C++ incorporee via Rcpp (sans parallelisation)
CalcMS4 <- function(imgAsMat){
  M <- foreach(i=1:3, .combine=c) %do% moy(imgAsMat[,i])
  #S <- foreach(i=1:3, .combine=c) %do% ect(imgAsMat[,i]-M[i]) # Ancienne version ect
  S <- foreach(i=1:3, .combine=c) %do% ect(imgAsMat[,i])
  return(cbind(M,S))
}

## Version C++ parallelisee dans R avec un mclapply
CalcMS5 <- function(imgAsLst) {
  M <- unlist(mclapply(imgAsLst, moy, mc.preschedule=FALSE, mc.cores=4))
  S <- unlist(mclapply(imgAsLst, ect, mc.preschedule=FALSE, mc.cores=4))
  return(cbind(M,S))
}

## Version C++ parallelisee dans C++ (via openMP)
CalcMS6 <- function(imgAsLst) {
  M <- unlist(lapply(imgAsLst, moy2))
  S <- unlist(lapply(imgAsLst, ect2))
  return(cbind(M,S))
}


#### Importing and executing PCA on JPEG files ####

#' PCA on the RGB channels of a JPEG file
#' 
#' JPEG images are converted into a n*3 matrix, where n stands for the number of pixels in the image and
#' the columns are the Red, Green and Blue (RGB) channels of the pixels. Then, a PCA is performed.
#' Outputs are two JPEG images obtained by projecting the original one to the (two) principal component(s).
#' @param file Character string specifying the path to the JPEG file to analyse
#' @return  Produces two JPEG files named img1.jpeg and img2.jpeg.
#' img1.jpeg is obtained from file by conserving only the first principal component.
#' img2.jpeg is obtained from file by conserving only the two first principal component(s).
#' @examples
#' data(lenna)
#' jpeg::writeJPEG(lenna, target='lenna.jpeg')
#' pcajpeg(file="lenna.jpeg")
#' @export
pcajpeg <- function(file = NULL){
  if (!is.character(file) & !is.null(file)) {stop("File must be a character string specifying the path to JPEG file")}
  else {}
  if (is.null(file)) {file = file.choose()} else {}
  img <- readJPEG(file)
  ## Taille de l'image
  d <- dim(img)
  n <- d[1] #largeur
  m <- d[2] #Hauteur
  ## Linearisation de l'image
  R <- as.numeric(img[,,1])
  G <- as.numeric(img[,,2])
  B <- as.numeric(img[,,3])
  ## Representation de l'image sous forme d'une matrice avec n*m lignes et 3 colonnes
  imgAsMat <- cbind(R, G, B)
  ##Calcul de la moyenne et de l'ecart-type de chaque canal
#   t <- proc.time()
  MS <- CalcMS2(imgAsMat)
#   print(proc.time()-t)
  ## Decommenter le prochain bloc d'instructions pour une optimisation du calcul de la moyenne et de l'ecart-type
#   L <- list(R,G,B)
#   t <- proc.time()
#   MS <- CalcMS6(L)  
#   print(proc.time()-t)
  print(MS)
  M <- MS[,1]
  S <- MS[,2]
  ## Realisation de l'ACP
  pca.img <- PCA(imgAsMat)
  # pca.img2 <- PCA(scale(imgAsMat), scale.unit = TRUE, graph=FALSE)
  # Decommenter la ligne suivante pour executer l'ACP sans representation graphique
#   pca.img <- PCA(imgAsMat, graph=FALSE)
  # Pass <- eigen(cor(imgAsMat))$vectors # Facteurs principaux
  Pass <- pca.img$var$coord
  Comp <- pca.img$ind$coord
  Comp <- scale(Comp)
  ## Reconstitution de l'image a partir de la premiere composante
  tmp <- cbind(Comp[,1], rep(0,n*m), rep(0, n*m))
  tmp <- tmp%*%t(Pass)
  tmp <- tmp*matrix(S, n*m, 3, byrow=TRUE) + matrix(M, n*m, 3, byrow=TRUE)
  img1 <- array(tmp, dim=d)
  ## Enregistrement de l'image dans le repertoire de travail
  writeJPEG(img1, target='img1.jpeg', quality=1)
  #### Reconstitution de l'image a partir des deux premieres composantes ####
  tmp <- cbind(Comp[,1:2], rep(0,n*m))
  tmp <- tmp%*%t(Pass)
  tmp <- tmp*matrix(S, n*m, 3, byrow=TRUE) + matrix(M, n*m, 3, byrow=TRUE)
  img2 <- array(tmp, dim=d)
  writeJPEG(img2, target='img2.jpeg', quality=1)
  #### Reconstitution de l'image a partir des trois composantes ####
  tmp <- pca.img$ind$coord
  tmp <- scale(tmp)%*%t(Pass)
  tmp <- tmp*matrix(S, n*m, 3, byrow=TRUE) + matrix(M, n*m, 3, byrow=TRUE)
  img3 <- array(tmp, dim=d)
  writeJPEG(img3, target='img3.jpeg', quality=1)
  print('Done!')
  #return(NULL)
}
#' @useDynLib acpjpeg
#' @import foreach
#' @importFrom Rcpp sourceCpp
NULL
#' @importFrom jpeg readJPEG
NULL
#' @importFrom jpeg writeJPEG
NULL
#' @importFrom FactoMineR PCA
NULL
#' @importFrom parallel mclapply
NULL

