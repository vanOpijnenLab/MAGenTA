#PLOT TRACKS
getTracksNew<-function(fromCoord,toCoord,x1,x2){
  #Fill in values if surrounding region coordinates weren't specified
  if (missing(x1))x1=0
  if (missing(x2))x2=0
  
  #HIGHLIGHT TRACK: encompasses several tracks
  ht1 <- HighlightTrack(trackList = list(dataTrack1,dataTrack2), start =fromCoord,end= toCoord,chromosome = "genome")
  
  plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkblue",fontsize=17,from=fromCoord-x1,to=toCoord+x2,main="Region of Interest")
  
}

getTracks<-function(fromCoord,toCoord,x,y){
  #Load libraries
  library(seqinr)
  library(Biostrings)
  library(Gviz)
  
  genomeFile="NC_012469.fasta"
  geneFile="19F_01genes.txt"
  singleFitFile="singleFit_19F_SDDM.txt"
  
  
  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack<-SequenceTrack(genomeFile,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF<-as.data.frame(read.table(file = geneFile,header=TRUE))
  anTrack<-AnnotationTrack(start=geneDF$start,end=geneDF$end,chromosome="genome",strand=geneDF$strand,id=geneDF$id,showFeatureId=TRUE, name="Genes",fill="gray",fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit<-as.data.frame(read.table(file=singleFitFile))
  insertFit$seqnames="genome"
  insertFit<-insertFit[c("seqnames","V1","V2","V3")]
  colnames(insertFit)<-c("seqnames","start","end","fit")
  insertFit$fit=insertFit$fit-1;
  insertFit$seqnames="genome"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,ignore.strand=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="genome",ylim=c(-1,1),type="histogram",
                      col.baseline="black",baseline=0,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(insertFit,fromCoord,toCoord)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText<-paste0("Window Average= ",avg,"  |   Actual Aggregate Fitness= ", avg+1)
  avgTrack<-AnnotationTrack(start=fromCoord,end=toCoord,chromosome="genome",id=avgText,showFeatureId=TRUE,fill="darkred", name="Window average")
  
  #Track for individual fitness values when the window is small enough
  fitTrack<-AnnotationTrack(insertFit,chromosome="genome",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
#Keep track pollution low:
  print(paste0("# of insertions in the window: ", count))
  if (count>15 || (toCoord-fromCoord>10000)){
    ht1 <- HighlightTrack(trackList = list(aggTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
  }
  else{
    ht1 <- HighlightTrack(trackList = list(aggTrack,fitTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
    
  }
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
             #background.title="darkred",fontsize=17,chromosome="genome",
             #from=fromCoord,to=toCoord)
  
  #TO CHANGE SURROUNDING WINDOWS ALTER:
  if (missing(x))x=50
  if (missing(y))y=50

  #plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y,main="19F from Galaxy")
  plotTracks(list(gTrack,sTrack,anTrack,fitTrack),chromosome="genome",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y)
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord," to ",toCoord," with a window (",fromCoord-x,",",toCoord+y,")")
  return(confirm)
  
  
}

Wavg<-function(fit,fromCoord,toCoord){
  sum=0; count=0;i=1;
  while (fit[i,2]<=fromCoord) i<-i+1;
  
  while (fit[i,2]<=toCoord){
    print(fit[i,])
    sum<-sum+fit[i,4]; count<-count+1;
    i<-i+1;
  }
  avg<-(sum/count);
  output<-list(avg,count)
  return(output);
}


getTracks2<-function(fromCoord,toCoord,x,y){
  #Load libraries
  library(seqinr)
  library(Biostrings)
  library(Gviz)
  
  genomeFile="9-daptomycin/0-genome/19F_012469.fasta"
  geneFile="9-daptomycin/0-genome/19F_genes.txt .txt"
  #singleFitFile="singleFit_tigr4_SDDM.txt"
  singleFitFile="9-daptomycin/singleFit/19FDapto_singleFit.txt "
  
  #genomeFile="NC_012469.fa"
  #geneFile="19F_genes.txt"
  #singleFitFile="singleFit_19F_SDDM.txt"
  
  
  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack<-SequenceTrack(genomeFile,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF<-as.data.frame(read.table(file = geneFile,header=TRUE))
  colnames(geneDF)<-c("seqnames","start","end","width","strand","id")
  anTrack<-AnnotationTrack(start=geneDF$start,end=geneDF$end,chromosome="genome",strand=geneDF$strand,id=geneDF$id,showFeatureId=TRUE, name="Genes",fill="gray",fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit<-as.data.frame(read.table(file=singleFitFile))
  insertFit$seqnames="genome"
  insertFit<-insertFit[c("seqnames","V1","V2","V3")]
  colnames(insertFit)<-c("seqnames","start","end","fit")
  insertFit$fit=insertFit$fit;
  insertFit$seqnames="genome"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,ignore.strand=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="genome",ylim=c(-.5,.5),type="histogram",
                      col.baseline="black",baseline=0,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(insertFit,fromCoord,toCoord)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText<-paste0("Window Average= ",avg,"  |   Actual Aggregate Fitness= ", avg)
  avgTrack<-AnnotationTrack(start=fromCoord,end=toCoord,chromosome="genome",id=avgText,showFeatureId=TRUE,fill="darkred", name="Window average")
  
  #Track for individual fitness values when the window is small enough
  fitTrack<-AnnotationTrack(insertFit,chromosome="genome",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
  #Keep track pollution low:
  #print(paste0("# of insertions in the window: ", count))
  if (count>15 || (toCoord-fromCoord>10000)){
    ht1 <- HighlightTrack(trackList = list(aggTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
  }
  else{
    ht1 <- HighlightTrack(trackList = list(aggTrack,fitTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
    
  }
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
  #background.title="darkred",fontsize=17,chromosome="genome",
  #from=fromCoord,to=toCoord)
  
  #TO CHANGE SURROUNDING WINDOWS ALTER:
  if (missing(x))x=50
  if (missing(y))y=50
  
  plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y,main="PyR Regulon SP in Penn")
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord," to ",toCoord," with a window (",fromCoord-x,",",toCoord+y,")")
  return(confirm)
  
}