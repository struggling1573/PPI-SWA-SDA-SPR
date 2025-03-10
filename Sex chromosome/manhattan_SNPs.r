library(qqman)
library(Cairo)

args=commandArgs(TRUE);
infile  <- args[1];
outfile  <- args[2];

Fstfile<-read.table(infile, header=F, stringsAsFactors=F)
SNP<-paste(Fstfile[,1],Fstfile[,2],sep = ":")
Fstfile=cbind(SNP,Fstfile)
colnames(Fstfile)<-c("SNP","CHR","POS","Fst")
Fstfile$SNP<-as.character(Fstfile$SNP)
Fstfile$CHR<-as.integer(Fstfile$CHR)
Fstfile<-na.omit(Fstfile)
filePNG<-paste(outfile,"manhattan.png",sep=".")
CairoPNG(file=filePNG,width=1500,height=500)
colorset<-c("#FF0000","#FFD700","#2E8B57","#7FFFAA","#6495ED","#0000FF","#FF00FF")
manhattan(Fstfile, chr="CHR", bp="POS", p="Fst", snp="SNP", col=colorset, logp=F, suggestiveline = F, genomewideline = F, ylab="Fst", ylim=c(0,1), font.lab=4,cex.lab=1.2, main="plot", cex=0.8)

### plot single chromosome
if (FALSE){
    manhattan(subset(Fstfile,CHR=="1"),chr="CHR",bp="POS",p="Fst",snp="SNP", col=colorset,logp=F,suggestiveline = F, genomewideline = F,ylab="Fst",ylim=c(0,1),font.lab=4,cex.lab=1.2,"chr1",cex=0.8)
}

dev.off()
