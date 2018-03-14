
setwd("C:/Users/hslansd/Desktop/")
#---Reading Necessary Files
allUTRs<-read.csv("All_UTRLengths.csv", header = TRUE, row.names = 1)
allmirs<-read.csv("All_Mirs.csv", header = TRUE, row.names = 1)

single_mirs<-read.csv("PainPTSDmirs.csv", header = TRUE, row.names = 1)

UTRs=allUTRs$transcript
#---Creating simulations
mc=replicate(n=10000,expr=sample(x=UTRs, size =402, replace = F))



#---Intializing a matrix to hold the occurences of each miR seen in each simulation (miRs not on pain list take out in Excel after R)

mc_matrix<-matrix(data=NA, nrow= 583, ncol=10000)
temp<-allmirs[allmirs$Transcript.ID %in% mc[,1],]
row.names(mc_matrix)<-temp_table$Var1

#---Creating MC_matrix

for(i in 1:10000){ 
  temp<-allmirs[allmirs$Transcript.ID %in% mc[,i],] #correct 
  row.names(temp)<-c(1:nrow(temp))
  temp_table<-as.data.frame(table(temp$miR.Family))
  mc_matrix[,i]<-temp_table$Freq
  }
write.csv(mc_matrix,"MC_MirMatrix.csv")

#---Collecting the UTR length for each gene in each simulation, 2615.45 is known from outside calculation of Pain/PTSD gene average length
UTRlength_matrix<-matrix(data=NA, nrow=402, ncol=10000)

for(i in 1:10000){
  temp<-allUTRs[allUTRs$transcript %in% mc[,i],]
  UTRlength_matrix[,i]<-temp$length}

#---Creating a scale factor to normalize miR count according to length
UTRaverages<-colSums(UTRlength_matrix)/402
UTRnorm<-2615.45/UTRaverages

#---Reforming matrix according to scale
norm_matrix<-matrix(data=NA,nrow=583,ncol=10000 )
  
for(i in 1:10000){norm_matrix[,i]<-mc_matrix[,i]*UTRnorm[i]}  
row.names(norm_matrix)<-row.names(mc_matrix)

write.csv(norm_matrix,"MCMatrix_NormCounts.csv")

#---Defining p values MCMatrix_NormCounts.csv has had non pain/ptsd mirs removed, and the occurencs of the 402 pain/ptsd genes placed in column 1

Full_Matrix<-read.csv("MCMatrix_NormCounts.csv", header=TRUE, row.names = 1)
f_mat<-as.matrix(Full_Matrix)

avgs<-rowMeans(f_mat)
sds<-apply(f_mat,1,sd)
z_scores<-as.numeric(list())
for(i in 1:560){
  z_scores[i]<-(f_mat[i,1]-avgs[i])/sds[i]}

p_values<-2*pnorm(-abs(z_scores))

