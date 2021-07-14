##################################################
#
# aim: pairwise FST by population table and plot starting from vcf file
# needs: "SNPRelate" Rpackage, vcf input file, population metadata tab separated input file
# notes: Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and are automatically replaced with zero.
#
# Author: Vincenzo Landi, University of Bari Aldo Moro, Italy
#vincenzo.landi@uniba.it
#modified from https://github.com/jiyongkun/pairwise_Fst/blob/master/Fst.R
# Date: 2021
##########################################################################################################################################
library("SNPRelate")
library(readr)

# need: 1) ped/map files or VCF ; 2) FAM or PED file for metafile preapration

##############editable part###############################################################
#----------conversion to gds file--------------------------------------------------------#
# edit only these lines (at line 138 you can change color of plot)#
ped.fn <- "./../patch_to_file_ped/fileped.ped" # pach to ped
map.fn <- "./../patch_to_file_map/filemap.map" #pach to map
#vcf.infile<-"patch_to_vcf/file_vcf.vcf" #¶pach to vcf
#--------------------retrieve FIF and IID from FAM or PED file
fam <- read_table2("./../patch_to_file_ped/fileped.fam",   col_names = FALSE)
#fam <- read_table2("./../patch_to_file_ped/fileped.ped",   col_names = FALSE)#




###########################################################################################
###########################################################################################



#-----------------------------------------------------------------------------------------#
snpgdsPED2GDS(ped.fn, map.fn, "ccm.gds") # from plink to gds
#snpgdsVCF2GDS(vcf.infile, "ccm.gds",  method="biallelic.only")# from vcf to gds

#-----------------load gds file-----------------------------------------------------------#

genofile.infile <- "ccm.gds"
transformNegativeFSTtozero<- TRUE
keepLowerTriangularMatrixOnly<- FALSE # change to TRUE to build only one half matrix
genofile <- snpgdsOpen(genofile.infile)

#---------------------metadata file ------------------------------------------------------#
metadata=fam[,c(2,1)]
colnames(metadata)=c("sample", "pop")



pop_code=metadata$pop
sample.id=read.gdsn(index.gdsn(genofile, "sample.id"))
poplevels=unique(pop_code)

#----------------pairwise FST calculation------------------------------------------------#


# pairwise populations matrix creation
res= outer(X= poplevels , Y= poplevels, 
           FUN = function(h,k){
             paste(h,k,sep = "/")
           }
)
colnames(res)=poplevels
rownames(res)=poplevels

as.data.frame(res)

for(i in poplevels) {
  for(j in poplevels) {
    if(grepl("/", res[i,j])) { #判断矩阵某一格是否包含"/"，从而只计算那些尚未计算的population pairs
      popelem= unlist(strsplit(res[i,j],"/"))
      
      #takes selection of samples an population to use for each pair
      flag<- pop_code %in% c(popelem[1],popelem[2])
      samp.sel<- sample.id[flag]
      pop.sel<- pop_code[flag]
      
      if (popelem[1]==popelem[2]){result="0"}else{
        results = snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), 
                            autosome.only=FALSE, method="W&C84",verbose=TRUE)
        result = results$Fst
      }
      res[i,j]=as.character(result)
      res[j,i]=as.character(result) #copy the result
    }
  }
}
#----------------final edits and prints------------------------------------------------#



# cell transformation from characters to numerics

res=data.frame(apply(res, 2, function(x) as.numeric(as.character(x))))
rownames(res)=poplevels

# trasforming negative FST into zero
# Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and will be automatically replaced with zero inside this function.
if (transformNegativeFSTtozero==TRUE){
  for(i in colnames(res)) {
    for(j in colnames(res)) {
      if (res[i,j]<0){res[i,j]<-0} 
    }
  }
}

# keep only the lower triangular matrix
# (set the upper triangular to zero)
if(keepLowerTriangularMatrixOnly == TRUE){
  res[upper.tri(res)]<-0
}
#-------------------------------------------------------------------------------------------#######
# PRINT pairwise FST matrix to text file
########
timewithspaces= Sys.Date()

timeAttr= gsub(pattern = "-",replacement = "_",x =  timewithspaces)

outfile <- paste("pairFstMatrix_", timeAttr, ".txt", sep="")

write.table(x = res, file = outfile, sep = "\t", dec = ".", quote = F, row.names = T,col.names = NA)


#---------------------edited from Arlequin's pairFstMatrix.r (Author: Heidi Lischer)--------#


numericMatrix=res

# preliminar functions
#----Mirror matrix (left-right)---------------------------------------------------------------#
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

#----Rotate matrix 270 clockworks-----------------------------------------------------------#
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

Matrix <- rotate270.matrix(numericMatrix)


# choose blue or orange gradient
ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))
#ColorRamp <- colorRampPalette(c("white",'orange1',"#8B4513"))


#------------------------- set name according to pc date---------------------------------#
timewithspaces= Sys.Date()
timeAttr= gsub(pattern = "-",replacement = "_",x =  timewithspaces)
outfileGraphic <- paste("pairFstMatrix_", timeAttr, ".png", sep="")
outfileGraphic_pdf <- paste("pairFstMatrix ", timeAttr, ".pdf", sep="")




#save graphic
#png(outfileGraphic, width=1300, height=1300, res=144)
pdf(outfileGraphic_pdf, width = 10, height = 10)  

smallplot <- c(0.874, 0.9, 0.18, 0.83)
bigplot <- c(0.13, 0.85, 0.14, 0.87)

old.par <- par(no.readonly = TRUE)

# draw legend --------------------------------
par(plt = smallplot)

# get legend values
Min <- min(Matrix, na.rm=TRUE)
Max <- max(Matrix, na.rm=TRUE)
binwidth <- (Max - Min) / 64
y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
z <- matrix(y, nrow = 1, ncol = length(y))

image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)

# adjust axis if only one value exists
if(Min == Max){
  axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
} else {
  axis(side=4, las = 2, cex.axis=0.8)
}

box()
mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)


#draw main graphic ---------------------------
a <- ncol(numericMatrix)
b <- nrow(numericMatrix)

x <- c(1:a)
y <- c(1:b)

par(new = TRUE, plt = bigplot)

image(x,y,as.matrix(Matrix), col=ColorRamp(64),
      main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",
      ylab="", axes=FALSE)
box()

#add labels
Labels=poplevels
if(is.null(Labels)){
  axis(1, at = c(1:a))
  axis(2, at = c(1:b), labels=c(b:1))
  mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
        font=2)
  mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
        font=2)
} else{
  axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
       las=2)
  axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
       las=2)
}


par(old.par)  #reset graphic parameters

dev.off()
############## END FSTPAIRWISE TABLE PLOT ############################################
#-------------------close genofile to clear enviroment-------------------------------#
snpgdsClose(genofile)


