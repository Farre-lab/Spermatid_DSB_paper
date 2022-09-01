#R code to describe how the stackplot of the coverage of states per chromosome was generated.

#Load the required libraries:
library(ggplot2)

#Set the working directory:
setwd("C:/Users/path/to/input/files")

#Read in the input file from ChromHMM- this file already contains the the longer stretches in the same state:
ST_16_states <- read.table("ST_16_segments.bed", header=FALSE)

#Make a length column:
ST_16_states$len <- ST_16_states$V3- ST_16_states$V2 #calculate length of each row

#Calculate the genome size:
genomeSize<-sum(ST_16_states$len) #I assume all genome has a state

coverage<-ST_16_states %>% group_by(V4) %>% summarise(percent=(sum(len)/genomeSize*100)) #calculate coverage of each state genomewide


chrSizeAll<-as.data.frame(ST_16_states %>% group_by(V1) %>% summarise(size=sum(len)))

#The input chromHMM file contains random/unplaced/chrM regions- this code makes sure they are removed from the data file
chrSizeAlmost<-chrSizeAll[!grepl("Un",chrSizeAll$V1),]
chrSize<-chrSizeAlmost[!grepl("random",chrSizeAlmost$V1),]
dataAll<-as.data.frame(ST_16_states %>% group_by(V1, V4) %>% summarise(totalLen=sum(len), .groups="keep"))
dataAlmost<-dataAll[!grepl("Un",dataAll$V1),]
data1<-dataAlmost[!grepl("M",dataAlmost$V1),]
data<-data1[!grepl("random",data1$V1),]

for(i in 1:nrow(chrSize)){
  print(chrSize[i,1])
  for(y in 1:nrow(data)){
    if(chrSize[i,1] == data[y,1]){
      data[y,4]<-data[y,3]/chrSize[i,2]*100
    }
  }
}

#Give each state a colour:
cols <- c("E1"="firebrick","E2"="coral","E3"="pale golden rod","E4"="yellow","E5"=" dark golden rod",
          "E6"="forest green", "E7"="medium aqua marine", "E8"="lime green", "E9"="light steel blue", "E10"="corn flower blue", "E11"="deep sky blue",
          "E12"="blue violet", "E13"="medium orchid", "E14"="purple", "E15"="plum", "E16"="deep pink")
          
 
#Create a factor for the chr & the states
X1 <- factor(data$V1, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13", "14","15","16","17","18", "19","20","21"))
V4.2 <- factor(data$V4, levels=c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13", "E14","E15","E16"))



#Plotting- from the input file I'd already changed the chr1 to 1 etc & chr X to 20 & chrY to 21.
#The X1 factor for the X axis makes the chr in order along the X-axis
#The V4.2 factor makes the states plot from E1 to E16 down the plot legend

pdf(" Stackplot ST 16 state chromHMM coverage per Chr.pdf", width = 10, height = 10, pointsize = 10)
ggplot(data, aes(fill=V4.2, y=V4.1, x=X1)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("% coverage of each chr by the 16 ST ChromHMM states") +
  scale_fill_manual(values=cols) +
  theme_classic() +  
  labs(x='Chr', y="% of the chr in each state", fill="ChromHMM state") 
dev.off() 

# On the X-axis 20 was then changed to X & 21 to Y in Inkscape.
