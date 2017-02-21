

  
# download needed packages

if(!require(cpgen)) {
  install.packages("cpgen"); require(cpgen)}

if(!require(rrBLUP)) {
  install.packages("rrBLUP"); require(rrBLUP)}

if(!require(matrixcalc)) {
  install.packages(matrixcalc); require(matrixcalc)}

if(!require(StAMPP)) {
  install.packages(StAMPP); require(StAMPP)}


#########Calculate the additive genomic relationship matrix

ARM<-function(M, type) {
     M2<-(M+1)/2  
     M3<-as.matrix(M)
  G<-switch(type,
         Yang = stamppGmatrix(M2),
         VanRaden = A.mat(M),
         Zhang = cgrm(M3)
  )
         rownames(G)<-rownames(M)
         colnames(G)<-rownames(M)
    return(G)
  }

  
#########Calculate the genomic relationship matrix inverse

G_inv<-function(GRM){
  GRM_pos<-pos.def(GRM)
  GRM_inv<-solve(GRM_pos)
}

##Test if positive semidefinite and fix


pos.def<-function(x) {
  if(is.positive.definite(x) == FALSE) { 
    make.positive.definite(x)
  } else
  { 
    return(x)
  }
}

#spectral decomposition to make positive semi-definite

make.positive.definite <- function(GRM)
{
  # M = GARM
  tol=1e-6
  eig <- eigen(GRM, symmetric=TRUE)
  rtol <- tol * eig$values[1]
  if(min(eig$values) < rtol)
  {
    vals <- eig$values
    vals[vals < rtol] <- rtol
    srev <- eig$vectors %*% (vals * t(eig$vectors))
    dimnames(srev) <- dimnames(GRM)
    return(srev)
  } else
  {
    return(GRM)
  }
}


#####Dominance Relationship Matrices

#Allele frequency table

allele.freq.func<-function(mm) { 
  markers.table<-lapply(mm, function(x) factor(x, levels=c(0,1,2)))
  alleles<-lapply(markers.table, table) 
  
  allele.freq.q<-sapply(alleles, function(x) (0.5*x[2]+x[1])/sum(x))
  allele.freq<-data.frame(allele.freq.q)
  colnames(allele.freq)<-"q"
  allele.freq$p<-1-allele.freq$q
  
  return(allele.freq)
}

######### dominance Centered relationship matrix

#calculate Dhat (centered): 

Dcent <- function(M) {
  M2<-M+1
  #make allele table
  allele_freq<-allele.freq.func(M2)
  
  #formulate the H matrix
  hcent<-function(g) { 
    
    p<-allele_freq$p
    q<-allele_freq$q
    
    hij <- ifelse(g == 0, -2*p*q,
                  ifelse(g == 1, 1 - 2*p*q, -2*p*q))
    return(hij)
  }
  h2<-apply(M2,2,hcent)
  
  #calculate the centering factor
  cent<-sum(apply(allele_freq,1, function(x) 2*x[1]*x[2]*(1-2*x[1]*x[2])))
  
  #get the cross product matrix
  drm_cent<-tcrossprod(h2)/cent
  #tcrossprod(replace(h2, is.na(h2), 0))/cent #from when missing values were replaced with a "zero",
  
  return(drm_cent)
}



######### dominance Normalized relationship matrix

#calculate Dhat (normalized): 

Dnorm <- function(M) {
  M2<-M+1
  #make allele table
  alleles<-allele.freq.func(M2)
  
    #formulate the H matrix
  hnorm<-function(g) { 
    
    p<-alleles$p
    q<-alleles$q
  
    hij <- ifelse(g == 0, -2*p^2,
                  ifelse(g == 0.5, 2*q*p, -2*q^2))
    return(hij)
  }
  
  h1<-apply(M2,2,hnorm)
  
  #calculate the normalizing factor
  norm<-sum(apply(alleles,1, function(x) (2*x[1]*x[2])^2))
  
  #get the cross product matrix
  tcrossprod(h1)/norm
  #tcrossprod(replace(h1, is.na(h1), 0))/norm
}

################

DRM<-function(M,type){
  switch(type,
    norm = Dnorm(M),
    center = Dcent(M)
  )  
}


########## TEST #############

#Create Data

lev<-c(-1,0,1) #assume biallelic markers with the notation "-1", "0" and "1"
markers <- data.frame(matrix(sample(lev,200*1000,replace=T),200,1000))

test1<-DRM(markers,"norm")
test2<-DRM(markers,"center")
test3<-ARM(markers,"Yang")
test4<-ARM(markers,"VanRaden")
test5<-ARM(markers,"Zhang")

