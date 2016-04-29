
#Create Data

lev<-c(0,0.5,1) #assume biallelic markers with the notation "0", "0.5" and "1"
M <- data.frame(matrix(sample(lev,200*1000,replace=T),200,1000))


############ Genomic Estimate Dominance Relationship

#Calculate the allele frequencies at each marker
# p is frequency for alleles "1" and half of the "0.5" alleles
# q is the frequency for the alleles "0" and half of the "0.5" alleles


markers.table<-lapply(M, function(x) factor(x, levels=c(0,0.5,1)))
alleles<-lapply(markers.table, table) 

allele.freq.q<-sapply(alleles, function(x) (0.5*x[2]+x[1])/sum(x))
allele.freq<-data.frame(allele.freq.q)
colnames(allele.freq)<-"q"
allele.freq$p<-1-allele.freq$q

############### Normalised Dominance Relationship Matrix
#calculate the H matrix (normalised) 
# This function is written to apply over columns:

hnorm<-function(g) { 
  
  #make the allele table
  markers.C<-factor(g,levels=c(0,0.5,1))
  alleles.C<-table(markers.C)
  allele.freq.q.C<-(0.5*alleles.C[2]+alleles.C[1])/sum(alleles.C)
  q<-(0.5*alleles.C[2]+alleles.C[1])/sum(alleles.C)
  p<-1-q
  
  #do the substitution
  hij <- ifelse(g == 0, -2*p^2,
                ifelse(g == 0.5, 2*q*p, -2*q^2))
  return(hij)
}

#apply the function across columns to the data:
h1<-apply(M,2,hnorm)


#Calculate the normalised dominance relationship:
# missing values are set as "0"
norm<-sum(apply(allele.freq,1, function(x) (2*x[1]*x[2])^2))

Dnorm <- function(x) {
  tcrossprod(replace(x, is.na(x), 0))/norm
}

#final genomic normalised dominance relationship matrix: 
DG.n<-Dnorm(h1)


########### Centered Dominance Relationship Matrix

# Calculate centered H matrix:

# This function is written to apply over columns:

hcent<-function(g) { 
  
  #make the allele table
  markers.C<-factor(g,levels=c(0,0.5,1))
  alleles.C<-table(markers.C)
  allele.freq.q.C<-(0.5*alleles.C[2]+alleles.C[1])/sum(alleles.C)
  q<-(0.5*alleles.C[2]+alleles.C[1])/sum(alleles.C)
  p<-1-q
  
  #do the substitution
  hij <- ifelse(g == 0, -2*p*q,
                ifelse(g == 0.5, 1 - 2*p*q, -2*p*q))
  return(hij)
}

#apply the function across columns to the data:
h2<-apply(M,2,hcent)


#calculate the centered dominance relationship matrix: 
# missing values are set as "0"

cent<-sum(apply(allele.freq,1, function(x) 2*x[1]*x[2]*(1-2*x[1]*x[2])))

Dcent <- function(x) {
  tcrossprod(replace(x, is.na(x), 0))/cent
}

#final genomic centered dominance relationship matrix: 
DG.c<-Dcent(h2)
