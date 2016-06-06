

############ Genomic Estimate Dominance Relationship

#Calculate the allele frequencies at each marker
# p is frequency for alleles "1" and half of the "0.5" alleles
# q is the frequency for the alleles "0" and half of the "0.5" alleles


#####Dominance Relationship Matrices

#Calcualte the Allele frequency table

allele.freq.func<-function(M) { 
  markers.table<-lapply(M, function(x) factor(x, levels=c(0,0.5,1)))
  alleles<-lapply(markers.table, table) 
  
  allele.freq.q<-sapply(alleles, function(x) (0.5*x[2]+x[1])/sum(x))
  allele.freq<-data.frame(allele.freq.q)
  colnames(allele.freq)<-"q"
  allele.freq$p<-1-allele.freq$q
  
  return(allele.freq)
}

######### Dominance Centered Relationship Matrix

#Calculate Dhat (centered): 

Dcent <- function(M) {
  
  #make allele table
  alleles<-allele.freq.func(M)
  
  #formulate the H matrix
  hcent<-function(g) { 
    
    p<-alleles$p
    q<-alleles$q
    
    hij <- ifelse(g == 0, -2*p*q,
                  ifelse(g == 0.5, 1 - 2*p*q, -2*p*q))
    return(hij)
  }
  
  h2<-apply(M,2,hcent)
  
  #calculate the centering factor
  cent<-sum(apply(alleles,1, function(x) 2*x[1]*x[2]*(1-2*x[1]*x[2])))
  
  #get the cross product matrix
  tcrossprod(replace(h2, is.na(h2), 0))/cent
}



######### dominance Normalized relationship matrix

#calculate Dhat (normalized): 

Dnorm <- function(M) {
  
  #make allele table
  alleles<-allele.freq.func(M)
  
    #formulate the H matrix
  hnorm<-function(g) { 
    
    p<-alleles$p
    q<-alleles$q
  
    hij <- ifelse(g == 0, -2*p^2,
                  ifelse(g == 0.5, 2*q*p, -2*q^2))
    return(hij)
  }
  
  h1<-apply(M,2,hnorm)
  
  #calculate the normalizing factor
  norm<-sum(apply(alleles,1, function(x) (2*x[1]*x[2])^2))
  
  #get the cross product matrix
  tcrossprod(replace(h1, is.na(h1), 0))/norm
}

################

#function to switch between either method of calculating Dhat

DRM<-function(M,type){
  switch(type,
    norm = Dnorm(M),
    center = Dcent(M)
  )  
}


########## TEST #############

#Create Data

lev<-c(0,0.5,1) #assume biallelic markers with the notation "0", "0.5" and "1"
markers <- data.frame(matrix(sample(lev,200*1000,replace=T),200,1000))

#test

test<-DRM(markers,"norm")
