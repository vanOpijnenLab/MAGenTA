#Scatter plot with conditional labeling

ggplot(dat, aes(x=mean.gluc, y=mean.dapto)) +
+     geom_point(shape=1,color="purple")+ geom_smooth(method=lm,se=TRUE) +geom_text(aes(label=ifelse(sub[,3]>.2,as.character(row.names(sub)),''),check_overlap=TRUE),size=3,hjust=-.2,vjust=-.4)+ggtitle("Gene fitness for TIGR4: Glucose v. Daptomycin")+xlab("Mean fitness for genes --glucose")+ylab("Mean fitness for genes --daptomycin")



###Whole process
table<-read.csv("19F_gluc-dapto.csv")
row.names(table)<-table[,1]
table[,1]<-NULL
x<-table[,1]
y<-table[,7]
dat<-data.frame(table)
genes<-row.names(dat)
dat[,3]<-abs(dat[,1]-dat[,2])


ggplot(dat, aes(x=mean.gluc, y=mean.dapto)) +geom_point(shape=1,color="red")+ geom_smooth(method=lm,se=TRUE) +geom_text(aes(label=ifelse(dat[,3]>.2,as.character(row.names(dat)),''),check_overlap=TRUE),size=3,hjust=-.2,vjust=-.4)+ggtitle("Gene fitness for 19F: Glucose v. Daptomycin")+xlab("Mean fitness for genes --glucose")+ylab("Mean fitness for genes --daptomycin")



ggplot(dat, aes(x=mean.gluc, y=mean.dapto)) +geom_point(shape=1,color="firebrick")+ geom_smooth(method=lm,se=TRUE,color="dodgerblue4") +geom_text(aes(label=ifelse(dat[,3]>.2,as.character(row.names(dat)),''),check_overlap=TRUE),size=3,hjust=1.2,vjust=-1)+ggtitle("Gene fitness for 19F: Glucose v. Daptomycin")+xlab("Mean fitness for genes --glucose")+ylab("Mean fitness for genes --daptomycin")+theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(face="bold", vjust=-0.35),axis.title.y = element_text(face="bold" , vjust=0.35))



"dodgerblue4", "darkolivegreen4",
"darkorchid3", "goldenrod1"
aquamarine4