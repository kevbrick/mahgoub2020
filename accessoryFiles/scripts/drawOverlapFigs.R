library(UpSetR)
library(nVennR)

source('genericFunctions.R')
theme7point()

x <- read.table('zcwpw1_overlaps.tab',header=TRUE)

usDF <- data.frame(nm=paste(x$cs,
                            x$from,
                            x$to,
                            sep=":"))
  
  
mylist <- vector(mode="list", length=4)

names(mylist) <- c('hs','tss','tes','cgi')
  
for (s in names(mylist)){
  mylist[[s]] <- usDF[x[[s]],1]
}

gStr <- ggplot(x,aes(x=hs,y=strength)) + 
  geom_violin(fill='grey90',lwd=.2) + 
  geom_boxplot(width=.05,outlier.size=.01,fill='pink',lwd=.2) + 
  xlab('Zcwpw1 peak overlaps hotspot') + ylab('Zcwpw1 CUT&RUN strength') + 
  coord_flip()

graphics.off()
png(filename = 'Zcwpw1_Strength_HS_v_others.png',
    width=3,height=3, res = 300, units='in')
print(gStr)
dev.off()

graphics.off()
svg(filename = 'Zcwpw1_Strength_HS_v_others.svg',
    width=3,height=3)
print(gStr)
dev.off()

## Show Ori Overlaps
graphics.off()
png(filename = 'Zcwpw1_UpSet.png',
    width=7,height=6, res = 300, units='in')
upset(fromList(mylist), order.by = "freq")
dev.off()
  
graphics.off()
svg(file = 'Zcwpw1_UpSet.svg',width=7,height=6)
upset(fromList(mylist), order.by = "freq")
dev.off()

vennSVG <- 'Zcwpw1_Venn.svg'
venn <- plotVenn(mylist[names(mylist)],
                 nCycles=5000,
                 systemShow=TRUE,
                 fontScale=1.25,
                 setColors=c('orange','grey','firebrick','dodgerblue'),
                 borderWidth=3)
showSVG(nVennObj = venn, outFile = vennSVG)
