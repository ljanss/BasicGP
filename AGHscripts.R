# Functions to compute A, G and H matrix as needed in BGLR (these are not the inverses!)
# Some examples in the bottom are commented, so they don't run when using source("AGHscripts.R")
# to load these functions.
# Needs to have package SparseM and pedigree installed.
# Usage:
# A <- makeAmat(ped)
#      where ped is a data table with 3 columns ID, parent1, parent2
# G <- makeGmat(geno, idlist, method=1)
#      where geno is table of genotypes as 0,1,2 (missings NA)
#      idlist is optional list of names for the IDs (will be appended on rows and columns when given)
#      method is standard set to 1, can also be set to 2 to get Van Raden 2 G-matrix
# H <- makeHmat(A, G, Aweight=0)
#      A: an A-matrix with names of individuals as the row-names of the matrix
#      G: a G-matrix with names of individuals as the row-names of the matrix
#         The names on the rows of A and G will be used to match genotyped and ungenotyped individuals
#      Aweight: optional weight (standard at 0) to weigh some part of A in an adjusted G matrix, when
#         used common values are 0.05 to 0.10.
#      The G-matrix is adjusted to the A-matrix using the Van Raden regression.

library(SparseM)
library(pedigree)

# A-matrix using package pedigree and package SparseM.
# Optional can indicate to make only part of A by setting which= with list of TRUE/FALSE for which IDs to make A,
# without setting of which= A is made for all IDs.
# Names (first column of ped) are put on rows and columns of A.
makeAmat <- function(ped, which=NULL, inbreeding=NULL) {
    if(is.null(which)) which=rep(TRUE,nrow(ped))
    if(is.null(inbreeding)) {
   	   makeA(ped,which = which)
   	}
   	else {
   	   makeA(ped,which = which,inbreeding=inbreeding)
   	}
	Alist <- scan("A.txt", what=list(integer(),integer(),double()),quiet=TRUE)
    A <- as.matrix(new("matrix.ssr",ra=Alist[[3]], ja=Alist[[2]], ia=c(which(!duplicated(Alist[[1]])),
                    (length(Alist[[1]])+1L)), dimension=c(sum(which),sum(which))))
    rownames(A) <- as.character(ped[which,1])
    colnames(A) <- as.character(ped[which,1])
    return(A)
}

# G-matrix with methods 1 and 2.
# If an idlist is given, names will be put on rows and columns.
makeGmat <- function(geno, idlist=NULL, method=1) {
	p <- colMeans(geno, na.rm=TRUE)/2
	twopqvec <- 2*p*(1-p)
	if(method==1) {
		geno <- scale(geno,scale=FALSE)
		scaling <- sum(twopqvec)
	}
	else if (method==2) {
		geno <- scale(geno,scale=sqrt(twopqvec))
		scaling <- ncol(geno)
	}
	else {
		cat(paste("Method not available: ",method,"\n",sep=""))
		return(NULL)
	}
	geno[is.na(geno)] <- 0
	G <- tcrossprod(geno) / scaling
	if(!is.null(idlist)) {
		rownames(G) <- as.character(idlist)
		colnames(G) <- as.character(idlist)
	}
	return(G)
}

# Make single-step H matrix from pedigree A and genomic G.
# This function assumes there are names on the rows of A and G to match and combine A and G.
# Also including a Van Raden style adjustment of G to A.
makeHmat <- function(A,G, Aweight=0, adjust=TRUE) {
   # indicate genotyped and non-genotyped (list of size A matrix)
   genotyped <- rep(FALSE,nrow(A))                     # first make a whole vector of FALSE
   geno2ped <- match(rownames(G),rownames(A))
   genotyped[geno2ped] <- TRUE                         # then where G matches A set TRUE
   cat(paste("Total in pedigree",nrow(A),"genotyped",sum(genotyped),"\n"))
   # parts and intermediary steps needed
   A12 <- A[!genotyped, genotyped]
   A22 <- A[genotyped, genotyped]
   A12A22inv <- A12 %*% solve(A22)
   # if G has unmatched IDs not present in A, reduce G to have only the matched ones
   if(sum(is.na(geno2ped))>0) {
      G <- G[!is.na(geno2ped),!is.na(geno2ped)]
   }
   # do Van Raden regression G ~ A22 to adjust G
   # this step went wrong when G has IDs that don't match A
   if (adjust) {
      fit <- lm(as.vector(G) ~ as.vector(A22))
      intercpt <- summary(fit)$coefficients[1,1]
      slope <- summary(fit)$coefficients[2,1]
      cat(paste("G to A22 adjustment has intercept",intercpt,"and slope", slope,"\n"))
      G <- (G-intercpt)/slope
   }
   # If Aweight >0 put some A-matrix in the G (this is from the A22 part)
   if (Aweight > 0) {
      G <- Aweight*A22 + (1-Aweight)*G
   }
   # make H matrix
   H <- matrix(0,nrow(A),nrow(A))
   H[!genotyped, !genotyped] <- A[!genotyped, !genotyped] + A12A22inv %*% (G - A22) %*% t(A12A22inv)
   H[!genotyped, genotyped] <- A12A22inv %*% G
   H[genotyped, !genotyped] <- t(H[!genotyped, genotyped])
   H[genotyped, genotyped] <- G
   rownames(H) <- rownames(A)
   colnames(H) <- rownames(A)
   return(H)
}

# example ped to make A
#id <- 1:6
#dam <- c(0,0,1,1,4,4)
#sire <- c(0,0,2,2,3,5)
#ped <- data.frame(id,dam,sire)
#makeAmat(ped)

# using the mouse pedigree
#mouseped <- read.table("mouseall_ped.txt")
#A <- makeAmat(mouseped[,1:3])

# making G matrix
#geno <- read.table("mousegeno_geno.txt")
#G <- makeGmat(geno[,-1],idlist=geno[,1])

# make H-matrix
#H <- makeHmat(A,G)


